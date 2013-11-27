/*!
 * 2D LU factorization
 *
 * Authors: Khalid Hasanov
 */

#include "lu_factorization.h"


inline double lu_factorize(LU_data* lu_data, Platform_data* platform_data) {

    double *B_a, *B_b; // matrix blocks
    double alpha = 1, beta = 1; // C := alpha * a * b + beta * c
    size_t B_proc_col, B_proc_row; // Number of bloc(row or col) on one processor

    int nb_proc;
    size_t size_row = platform_data->size_row;
    size_t size_col = platform_data->size_col;
    size_t k = lu_data->k_global;
    size_t k_a = k / size_row;
    size_t k_b = k / size_col;
    size_t m = lu_data->m_global / size_col;
    size_t n = lu_data->n_global / size_row;
    lu_data->m = m;
    lu_data->n = n;
    lu_data->k_a = k_a;
    lu_data->k_b = k_b;

    size_t Block_size = lu_data->Block_size;
    size_t key = lu_data->key;

    B_proc_col = k_b / Block_size; // Number of block on one processor
    B_proc_row = k_a / Block_size; // Number of block on one processor

    size_t lda = k_a;
    size_t ldb = n;
    size_t ldc = n;
    size_t lda_local = lda;
    size_t ldb_local = ldb;

    lu_data->nb_block = k / Block_size;

    init_timer();

    struct timespec start_time, end_time; //time measure
    struct timespec start_time_intern, end_time_intern; //time measure
    double time, communication_time = 0, computation_time;


    /*--------------------Communication types for MPI--------------------------*/
    MPI_Datatype Block_a;
    MPI_Datatype Block_a_local;
    MPI_Datatype Block_b;
    MPI_Type_vector(m, Block_size, k_a, MPI_DOUBLE, &Block_a);
    MPI_Type_vector(m, Block_size, Block_size, MPI_DOUBLE, &Block_a_local);
    MPI_Type_vector(Block_size, n, n, MPI_DOUBLE, &Block_b);
    MPI_Type_commit(&Block_a);
    MPI_Type_commit(&Block_a_local);
    MPI_Type_commit(&Block_b);
    /*-------------Communication types for MPI are configured------------------*/


    //  MPI_Barrier(platform_data->comm);

    MPI_Comm my_world;
    MPI_Comm_split(platform_data->comm, 1, key, &my_world);

    MPI_Comm_size(my_world, &nb_proc);
    MPI_Comm_rank(my_world, &platform_data->my_rank);


    if (nb_proc < (int) (size_row * size_col)) {
        info_print(0, (size_t) 0, (size_t) 0,
                "Not enough processors nb_proc : %d required : %zu\n",
                nb_proc, size_row * size_col);
        return -1;
    }


    size_t grid_proc = size_row*size_col;
    size_t row = platform_data->my_rank / size_row;
    size_t col = platform_data->my_rank % size_row;


    printf("my_rank:%d, size_row:%zu, size_col:%zu, row:%zu, col%zu\n",
            platform_data->my_rank, size_row, size_col, row, col);


    double *a, *b, *c;
    matrices_initialisation(&a, &b, &c, m, k_a, k_b, n, row, col);

    double *a_local, *b_local;
    blocks_initialisation(&a_local, &b_local, m, Block_size, n);


    if (lu_factorize(lu_data, platform_data) == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    if (platform_data->useless == 1) {
        /*----------------------Prepare the Communication Layer-------------------*/
        /* add useless processor on a new color to execute the matrix
         * multiplication with the other processors*/

        /* Split comm size_to row and column comms */
        MPI_Comm row_comm, col_comm;
        /* color by row, rank by column */
        MPI_Comm_split(my_world, size_row, MPI_UNDEFINED, &row_comm);
        /* color by column, rank by row */
        MPI_Comm_split(my_world, size_col, MPI_UNDEFINED, &col_comm);
        /*------------------------Communication Layer can be used-----------------*/

        return 0;
    }

    debug_print(1, row, col, "I'm initialized col: %zu row: %zu, "
            "size_col: %zu , size_row: %zu, my rank: %d \n"
            , col, row, size_col, size_row, platform_data->my_rank);



    /*------------------------Prepare the Communication Layer-------------------*/
    /* Split comm size_to row and column comms */
    MPI_Comm row_comm, col_comm;
    /* color by row, rank by column */
    MPI_Comm_split(my_world, row, col, &row_comm);
    /* color by column, rank by row */
    MPI_Comm_split(my_world, col, row, &col_comm);
    /*-------------------------Communication Layer can be used------------------*/



    get_time(&start_time);

    int start = 0;
    int end = lu_data->nb_block;

    /*-------------Distributed LU factorization algorithm-----------------*/
    size_t iter;
    for (iter = start; iter < end; iter++) {
        size_t pivot_row, pivot_col, pos_a, pos_b;

        // pivot row on processor layer
        pivot_row = (iter % size_col);
        pivot_col = (iter % size_row);

        //position of the block
        pos_a = (size_t) (iter / size_row) * Block_size;
        pos_b = (size_t) (iter / size_col) * ldb * Block_size;

#if DEBUG
        MPI_Barrier(row_comm);
        MPI_Barrier(col_comm);
#endif

        get_time(&start_time_intern);

        //Broadcast the row
        if (size_row > 1) {
            MPI_Datatype *Block;
            if (pivot_col != col) {
                B_a = a_local;
                lda_local = Block_size;
                debug_print(1, row, col, "recieve B_a %zu,%zu \n", m, Block_size);
                Block = &Block_a_local;
            } else {
                B_a = a + pos_a;
                lda_local = lda;
                debug_print(1, row, col, "sent B_a %zu,%zu \n", m, Block_size);
                Block = &Block_a;
            }

            hpnla_bcast(B_a, 1, *Block, pivot_col, row_comm, lu_data->bcast_algorithm);

        } else {
            B_a = a + pos_a;
            debug_print(1, row, col, "position of B_a: %zu \n", pos_a);
        }

        //Broadcast the col
        if (size_col > 1) {
            if (pivot_row == row) {
                B_b = b + pos_b;
                debug_print(1, row, col, "sent B_b Block_size: %zu, pos:%zu \n", ldb, pos_b);
            } else {
                B_b = b_local;
                debug_print(1, row, col, "recieve B_b %zu,%zu \n", Block_size, n);
            }

            hpnla_bcast(B_b, 1, Block_b, pivot_row, col_comm, lu_data->bcast_algorithm);

        } else {
            B_b = b + pos_b;
            debug_print(1, row, col, "position of B_b: %zu \n", pos_b);
        }

        get_time(&end_time_intern);
        communication_time += get_timediff(&start_time_intern, &end_time_intern);

#if DEBUG
        MPI_Barrier(row_comm);
        MPI_Barrier(col_comm);
#endif
        get_time(&start_time_intern);

        {

            SMPI_SAMPLE_GLOBAL(30, 0.1) {
                hpnla_gemm_execute(
                        lu_data->gemm, //user parameters: here equal NULL
                        CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n,
                        Block_size, alpha, B_a, lda_local, B_b, ldb_local, beta,
                        c, ldc);
            }
        }

        get_time(&end_time_intern);
        computation_time += get_timediff(&start_time_intern, &end_time_intern);

    }

#if DEBUG
    MPI_Barrier(row_comm);
    MPI_Barrier(col_comm);
#endif

    get_time(&end_time);
    time = get_timediff(&start_time, &end_time);

    check_result(c, a, b, m, n, k_a, k_b, row, col, size_row, size_col);

    printf("communication time: %le microseconds, " "computation time: %le microseconds\n", communication_time, computation_time);


    // close properly the pragram

    MPI_Type_free(&Block_a);
    MPI_Type_free(&Block_a_local);
    MPI_Type_free(&Block_b);

    SMPI_SHARED_FREE(a_local);
    SMPI_SHARED_FREE(b_local);


    SMPI_SHARED_FREE(a);
    SMPI_SHARED_FREE(b);
    SMPI_SHARED_FREE(c);

    MPI_Comm_free(&my_world);
    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);

    return time;
}

int validate_input(LU_data* lu_data, Platform_data* platform_data) {

    if (lu_data->k_global % lu_data->Block_size != 0) {
        info_print(0, lu_data->row, lu_data->col,
                "The matrix size has to be proportional to the number\
                of blocks: %zu\n", lu_data->nb_block);
        return EXIT_FAILURE;
    }

    if (platform_data->size_row > lu_data->nb_block || platform_data->size_col > lu_data->nb_block) {
        info_print(0, lu_data->row, lu_data->col,
                "Number of blocks is too small compare to the number of"
                " processors (%zu,%zu) in a row or a col (%zu)\n",
                platform_data->size_col, platform_data->size_row, lu_data->nb_block);
        return EXIT_FAILURE;
    }

    if (lu_data->nb_block % platform_data->size_row != 0 || lu_data->nb_block % platform_data->size_col != 0) {
        info_print(0, lu_data->row, lu_data->col, "The number of Block by processor is not an %s",
                "integer\n");
        return EXIT_FAILURE;
    }

    platform_data->useless = 0;
    if (lu_data->row >= platform_data->size_col || lu_data->col >= platform_data->size_row) {
        info_print(0, lu_data->row, lu_data->col, "I'm useless bye!!! col: %zu row: %zu, "
                "size_col: %zu , size_row: %zu \n",
                lu_data->col, lu_data->row, platform_data->size_col, platform_data->size_row);
        platform_data->useless = 1;

        return EXIT_SUCCESS;
    }


}