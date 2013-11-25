/*!
 * Block-cylic distributed LU factorization 
 * algorithm with MPI
 *
 * Authors: Khalid Hasanov
 *
 */

#include "tools/hpnla_debug.h"
#include "tools/hpnla_timer.h"
#include "communication/hpnla_bcast.h"
#include "Matrix_init.h"
#include "lu_factorization.h"
#include <stdlib.h>
#ifdef HDNLA_SMPI
#include <smpi.h>
#else
#define SMPI_SHARED_FREE free
#endif

double lu_factorize(LU_data* lu_data, Platform_data* platform_data) {
    double *a, *b, *c;
    int useless = 0;
    size_t group_row = lu_data->group / platform_data->size_group_row;
    size_t group_col = lu_data->group % platform_data->size_group_row;
    
    lu_data->group_row = group_row;
    lu_data->group_col = group_col;

    size_t pivot_group_col = 0;
    size_t pivot_group_row = 0;
    size_t pivot_col = 0;
    size_t pivot_row = 0;

    size_t start = 0;
    size_t end = 0;
    size_t pos_a = 0;
    size_t pos_b = 0;

    size_t B_proc_col, B_proc_row;      // Number of bloc(row or col) on one processor
    size_t B_group_col, B_group_row;    // Number of bloc(row or col) on one processor

    init_timer();

    double time, communication_time = 0;
    struct timespec start_time, end_time;               
    struct timespec start_time_intern, end_time_intern; 
    
    MPI_Comm my_world;

    if (lu_data->group >= platform_data->size_group_col * platform_data->size_group_row) {
        info_print(0, (size_t) 0, (size_t) 0,
                "Not enough group NB_groups : %zu my group id : %zu\n",
                platform_data->size_group_col * platform_data->size_group_row, lu_data->group);
        MPI_Comm_split(platform_data->comm, 0, lu_data->key, &my_world);
        return -1;
    } else {
        MPI_Comm_split(platform_data->comm, 1, lu_data->key, &my_world);
    }

    MPI_Comm_size(my_world, &platform_data->nb_proc);

    if (platform_data->nb_proc < (int) (platform_data->nb_requested_proc)) {
        info_print(0, (size_t) 0, (size_t) 0,
                "Not enough processors nb_proc : %d required : %zu\n", platform_data->nb_proc,
                platform_data->nb_requested_proc);
        return -1;
    }

        
    MPI_Comm group_comm, group_comm_tmp;
    MPI_Comm_split(my_world, lu_data->group, lu_data->key, &group_comm_tmp);

    MPI_Comm_rank(group_comm_tmp, &platform_data->my_rank);

    size_t row = platform_data->my_rank / platform_data->size_row;
    size_t col = platform_data->my_rank % platform_data->size_row;
    
    lu_data->row = row;
    lu_data->col = col;
    

    /* 
     * matrix sizes on one processor after distribution
     */
    lu_data->m = lu_data->m_global / (platform_data->size_col * platform_data->size_group_col);
    lu_data->n = lu_data->n_global / (platform_data->size_row * platform_data->size_group_row);
    lu_data->k_a = lu_data->k_global / (platform_data->size_row * platform_data->size_group_row);
    lu_data->k_b = lu_data->k_global / (platform_data->size_col * platform_data->size_group_col);

    size_t m = lu_data->m;
    size_t n = lu_data->n;
    size_t k_a = lu_data->k_a;
    size_t k_b = lu_data->k_b;

    if (row >= platform_data->size_col || col >= platform_data->size_row) {
        info_print(0, row, col,
                "I'm useless bye!!! col: %zu row: %zu, " "size_col: %zu , size_row: %zu \n",
                col, row, platform_data->size_col, platform_data->size_row);
        useless = 1;
    }

    if (useless == 1) {
        /*----------------------Prepare the Communication Layer------------------*/
        /* add useless processor on a new color to execute the matrix
         * multiplication with the other processors*/

        /* Split comm size_to row and column comms */
        MPI_Comm my_last_world;
        MPI_Comm_split(my_world, MPI_UNDEFINED, 0, &my_last_world);
        MPI_Comm_free(&my_world);

        return 0;
    }

   

    /*------------------------Prepare the Communication Layer------------------*/
    /* Split comm size_to row and column comms */
    MPI_Comm row_comm, col_comm, group_col_comm, group_row_comm;
    MPI_Comm my_last_world;
    MPI_Comm_split(my_world, 0, lu_data->key, &my_last_world);
    MPI_Comm_free(&my_world);

    MPI_Comm_split(my_last_world, lu_data->group, platform_data->my_rank, &group_comm);

    MPI_Comm_split(my_last_world,
            col + row * platform_data->size_row + platform_data->size_col * platform_data->size_row * group_row, group_col,
            &group_row_comm);
    MPI_Comm_split(my_last_world,
            col + row * platform_data->size_row + platform_data->size_col * platform_data->size_row * group_col, group_row,
            &group_col_comm);
    /* color by row, rank by column */
    MPI_Comm_split(group_comm, row, col, &row_comm);
    /* color by column, rank by row */
    MPI_Comm_split(group_comm, col, row, &col_comm);
    /*-------------------------Communication Layer can be used------------------*/

    size_t ldb = n;
    size_t lda = k_a;
    MPI_Datatype Block_a;
    MPI_Datatype Block_a_local;
    MPI_Datatype Block_b;


    MPI_Type_vector(m, lu_data->Block_size_out / platform_data->size_row, lda, MPI_DOUBLE,
            &Block_a);

    MPI_Type_vector(m, lu_data->Block_size_out / platform_data->size_row, lu_data->Block_size_out,
            MPI_DOUBLE, &Block_a_local);

    MPI_Type_vector(lu_data->Block_size_out / platform_data->size_col, n, ldb, MPI_DOUBLE,
            &Block_b);

    MPI_Type_commit(&Block_a);
    MPI_Type_commit(&Block_a_local);
    MPI_Type_commit(&Block_b);

    double * B_a;
    double * B_b;
    double * a_local;
    double * b_local;

    matrices_initialisation(&a, &b, &c, m, k_a, k_b, n,
            row + group_row * platform_data->size_col, col + group_col * platform_data->size_row);
    blocks_initialisation(&a_local, &b_local, m, lu_data->Block_size_out, n);

    B_proc_col = k_b / lu_data->Block_size_out; // Number of block on one processor
    B_proc_row = k_a / lu_data->Block_size_out; // Number of block on one processor
    B_group_col = B_proc_col * platform_data->size_col; // Number of block on one group
    B_group_row = B_proc_row * platform_data->size_row; // Number of block on one group
    
    /* for each group_block broadcast it between the groups and compute it */
    MPI_Barrier(my_last_world);
    get_time(&start_time);

    /*Communicate the block between the groups */
    int iter;
    int NB_block_group = lu_data->k_global / lu_data->Block_size_out;
    
    if (validate_input(lu_data, platform_data) == -1) {
        printf("validate failed\n");
        return -1;
    }
    
    /*--------------------Allocation of matrices block inside a group----------*/
  /*  double *a_in, *b_in;
    blocks_initialisation(&a_in, &b_in, m, lu_data->Block_size_in, n);
   * 
   * */

    /*--------------------Communication types for MPI--------------------------*/
 /*   MPI_Datatype *   Block_Summa_a;
    MPI_Datatype Block_Summa_a_group_local;
    MPI_Datatype Block_Summa_a_group_global;
    MPI_Datatype Block_Summa_a_local;
    MPI_Datatype Block_Summa_b;

    MPI_Type_vector(m, lu_data->Block_size_in, lu_data->Block_size_out, MPI_DOUBLE,
            &Block_Summa_a_group_local);
    MPI_Type_vector(m, lu_data->Block_size_in, k_a, MPI_DOUBLE,
            &Block_Summa_a_group_global);
    MPI_Type_vector(m, lu_data->Block_size_in, lu_data->Block_size_in, MPI_DOUBLE,
            &Block_Summa_a_local);
    MPI_Type_vector(lu_data->Block_size_in, n, n, MPI_DOUBLE, &Block_Summa_b);

    MPI_Type_commit(&Block_Summa_a_group_local);
    MPI_Type_commit(&Block_Summa_a_group_global);
    MPI_Type_commit(&Block_Summa_a_local);
    MPI_Type_commit(&Block_Summa_b);
  * 
  * */
    /*-------------Communication types for MPI are configured-------------------*/

    for (iter = 0; iter < NB_block_group; iter++) {
        // pivot on group layer

        debug_print(0, row + group_row * platform_data->size_col,
                col + group_col * platform_data->size_row,
                "iter : %d, B_group_row %zu, B_group_col %zu, " "B_proc_row %zu, B_proc_col %zu \n",
                iter, B_group_row, B_group_col, B_proc_row, B_proc_col);
        pivot_group_col = (iter / B_group_row);
        pivot_group_row = (iter / B_group_col);

        // Which data to send????
        // Which processor send?
        // All number of block between group should be proportional to the
        // number of processor on one row or col of a group
        int size_block_a = lu_data->Block_size_out / platform_data->size_row;
        int size_block_b = lu_data->Block_size_out / platform_data->size_col;

       
        //position of the block
        start = iter * lu_data->Block_size_out / lu_data->Block_size_in;
        end = (iter + 1) * lu_data->Block_size_out / lu_data->Block_size_in;
        pos_a = (iter % B_group_row) * size_block_a;
        pos_b = (iter % B_group_col) * size_block_b * ldb;

       
        get_time(&start_time_intern);
        //Broadcast the row
        size_t lda_local;

        if (platform_data->size_group_row > 1) {
            MPI_Datatype * Block;

            if (pivot_group_col != group_col) {
                B_a = a_local;
                lda_local = lu_data->Block_size_out;
                debug_print(0, row + group_row * platform_data->size_col,
                        col + group_col * platform_data->size_row,
                        "recieve group B_a col %zu pivot_col %zu\n",
                        col, col );
                Block = &Block_a_local;
            //    Block_Summa_a = &Block_Summa_a_group_local;
            } else {
                B_a = a + pos_a;
                lda_local = lda;
                debug_print(0, row + group_row * platform_data->size_col,
                        col + group_col * platform_data->size_row,
                        "sent group B_a col %zu pivot_col %zu\n",
                        col, col
                        );

                Block = &Block_a;
            }

            hpnla_bcast(B_a, 1, *Block, pivot_group_col, group_row_comm,
                    lu_data->bcast_algorithm);

        } else {
            B_a = a + pos_a;
            lda_local = lda;    
        }

        //Broadcast the col
        if (platform_data->size_group_col > 1) {
            if (pivot_group_row != group_row) {
                B_b = b_local;
            } else {
                B_b = b + pos_b;
            }

            hpnla_bcast(B_b, 1, Block_b, pivot_group_row, group_col_comm,
                    lu_data->bcast_algorithm);


        } else {
            B_b = b + pos_b;
            debug_print(0, row + group_row * platform_data->size_col,
                    col + group_col * platform_data->size_row,
                    "No group Bcast position of B_b size_group_col: %zu \n",
                    platform_data->size_group_col);
        }

        MPI_Barrier(my_last_world);

        get_time(&end_time_intern);
        communication_time += get_timediff(&start_time_intern, &end_time_intern);
    }
    
   printf("Until this point ok : %d\n", platform_data->my_rank);

    
    /* the work is finished*/
    MPI_Barrier(my_last_world);

    get_time(&end_time);
    time = get_timediff(&start_time, &end_time);

    MPI_Barrier(my_last_world);

    // close the resources
  /*  MPI_Type_free(&Block_Summa_a_group_local);
    MPI_Type_free(&Block_Summa_a_group_global);
    MPI_Type_free(&Block_Summa_a_local);
    MPI_Type_free(&Block_Summa_b);
   */
    MPI_Type_free(&Block_b);
    MPI_Type_free(&Block_a_local);
    MPI_Type_free(&Block_a);
    
  //  SMPI_SHARED_FREE(a_Summa);
  //  SMPI_SHARED_FREE(b_Summa);

    SMPI_SHARED_FREE(a_local);
    SMPI_SHARED_FREE(b_local);

    SMPI_SHARED_FREE(a);
    SMPI_SHARED_FREE(b);
    SMPI_SHARED_FREE(c);

    MPI_Comm_free(&group_comm_tmp);
    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&my_last_world);
    MPI_Comm_free(&group_comm);
    MPI_Comm_free(&group_row_comm);
    MPI_Comm_free(&group_col_comm);

    return EXIT_SUCCESS;
}

int validate_input(LU_data* lu_data, Platform_data* platform_data) {
    size_t k_a = lu_data->k_a;
    size_t k_b = lu_data->k_b;
    size_t row = lu_data->row;
    size_t col = lu_data->col;
    size_t group_row = lu_data->group_row;
    size_t group_col = lu_data->group_col;
    
   
    size_t NB_Block = lu_data->k_global / lu_data->Block_size_in;
        
    
    if (k_a % lu_data->Block_size_in != 0 || k_a % lu_data->Block_size_out) {
        info_print(0, row, col,
                "The matrix size has to be proportionnal to the number\
                of blocks: %zu\n",
                NB_Block);
        return -1;
    }

    if ((k_a < k_b ? k_a : k_b) < lu_data->Block_size_out
            || (k_a < k_b ? k_a : k_b) < lu_data->Block_size_in) {
        info_print(0, row, col,
                "K: %zu should be bigger than\
                Block_size_group : %zu and Block_size %zu\n",
                k_a < k_b ? k_a : k_b, lu_data->Block_size_out, lu_data->Block_size_in);
        return -1;
    }

    if (lu_data->Block_size_out % lu_data->Block_size_in != 0) {
        info_print(0, row, col,
                "The number of Block_size_group %zu should be\
                proportionnal to Block_size : %zu\n",
                lu_data->Block_size_out, lu_data->Block_size_in);
        return -1;
    }

    if ((lu_data->Block_size_out / lu_data->Block_size_in) % platform_data->size_row != 0
            || (lu_data->Block_size_out / lu_data->Block_size_in) % platform_data->size_col) {
        info_print(0, row, col,
                "The Number of block per group block : %d\
                should be proportional to size_row : %zu and size_col %zu\n",
                (int) (lu_data->Block_size_out / lu_data->Block_size_in), platform_data->size_row, platform_data->size_col);

        return -1;
    }


    if (platform_data->size_row > NB_Block || platform_data->size_col > NB_Block) {
        info_print(0, row, col,
                "Number of blocks is too small compare to the number of" "processors (%zu,%zu) in a row or a col (%zu)\n",
                platform_data->size_col, platform_data->size_row, NB_Block);

        return -1;
    }

    if (NB_Block % platform_data->size_row != 0 || NB_Block % platform_data->size_col != 0) {
        info_print(0, row + group_row * platform_data->size_col, col + group_col * platform_data->size_row,
                "The number of Block by processor is not an integer\n");
        return -1;
    }

}