/*!
 * Hierarchical Block Matrix Multiplication
 *
 * Authors: Quintin Jean-Noël
 *
 */

#include "tools/hpnla_debug.h"
#include "himpi/mpi_bcast_algs.h"
#include "Matrix_init.h"
#include "Summa.h"
#include "tools/hpnla_timer.h"
#include <stdlib.h>
#ifdef HMPI_SMPI
#include <smpi.h>
#else
#define SMPI_SHARED_FREE free
#endif

double recursive(hpnla_gemm* gemm, size_t m, size_t k, size_t n,
		size_t Block_size, size_t Block_size_group, size_t group, size_t key,
		size_t size_row, size_t size_col, size_t size_group_row,
		size_t size_group_col, MPI_Comm world, int bcast_algorithm,
		int distribution) {
	double *a, *b, *c;
	double summa_computation_t = 0, summa_communication_t = 0;
	/* Split the communicator into groups */

	int is_cyclic = (distribution > 0);

	/* Find out my identity in the default communicator */
	int myrank;
	int NB_proc;
//  int err;
	int useless = 0;
	size_t global_m = m;
	size_t global_n = n;
	size_t global_k = k;
	size_t group_row = group / size_group_row;
	size_t group_col = group % size_group_row;

	size_t pivot_group_col = 0;
	size_t pivot_group_row = 0;
	size_t pivot_col = 0, pivot_row = 0;

	size_t start = 0;
	size_t end = 0;
	size_t pos_a = 0;
	size_t pos_b = 0;

	size_t B_proc_col, B_proc_row; // Number of bloc(row or col) on one processor
	size_t B_group_col, B_group_row; // Number of bloc(row or col) on one processor

	init_timer();

	double time, communication_time = 0;
	struct timespec start_time, end_time; 				//time mesure
	struct timespec start_time_intern, end_time_intern; //time mesure

	MPI_Comm my_world;

	if (group >= size_group_col * size_group_row) {
		info_print(0, (size_t ) 0, (size_t ) 0,
				"Not enough group NB_groups : %zu my group id : %zu\n",
				size_group_col * size_group_row, group);
		MPI_Comm_split(world, 0, key, &my_world);
		return -1;
	} else {
		MPI_Comm_split(world, 1, key, &my_world);
	}

	MPI_Comm_size(my_world, &NB_proc);

	if (NB_proc
			< (int) (size_row * size_col * size_group_col * size_group_row)) {
		info_print(0, (size_t ) 0, (size_t ) 0,
				"Not enough processors NB_proc : %d required : %zu\n", NB_proc,
				size_row * size_col * size_group_col * size_group_row);
		return -1;
	}

	MPI_Comm group_comm, group_comm_tmp;
	MPI_Comm_split(my_world, group, key, &group_comm_tmp);

	MPI_Comm_rank(group_comm_tmp, &myrank);
	MPI_Comm_size(group_comm_tmp, &NB_proc);
	/* for each group start the execution of his */

	NB_proc = size_row * size_col;
	size_t row = myrank / size_row;
	size_t col = myrank % size_row;

        printf("my_rank=%d, size_row=%zu, row=%zu, col=%zu\n", myrank, size_row, row, col);

        
        
	// matrix sizes
	m = m / (size_col * size_group_col);
	n = n / (size_row * size_group_row);
	size_t k_a = k / (size_row * size_group_row);
	size_t k_b = k / (size_col * size_group_col);

	/*-------------------------Check some mandatory conditions-----------------*/
	size_t NB_Block = k / Block_size;
	if (k_a % Block_size != 0 || k_a % Block_size_group) {
		info_print(0, row, col,
				"The matrix size has to be proportionnal to the number\
                of blocks: %zu\n",
				NB_Block);
		return -1;
	}

	if ((k_a < k_b ? k_a : k_b) < Block_size_group
			|| (k_a < k_b ? k_a : k_b) < Block_size) {
		info_print(0, row, col,
				"K: %zu should be bigger than\
                Block_size_group : %zu and Block_size %zu\n",
				k_a < k_b ? k_a : k_b, Block_size_group, Block_size);
		return -1;
	}

	if (Block_size_group % Block_size != 0) {
		info_print(0, row, col,
				"The number of Block_size_group %zu should be\
                proportionnal to Block_size : %zu\n",
				Block_size_group, Block_size);
		return -1;
	}

	if (is_cyclic) {
		if ((Block_size_group / Block_size) % size_row != 0
				|| (Block_size_group / Block_size) % size_col) {
			info_print(0, row, col,
					"The Number of block per group block : %d\
                should be proportional to size_row : %zu and size_col %zu\n",
					(int )(Block_size_group / Block_size), size_row, size_col);

			return -1;
		}
	}

	if (size_row > NB_Block || size_col > NB_Block) {
		info_print(0, row, col,
				"Number of blocks is too small compare to the number of" "processors (%zu,%zu) in a row or a col (%zu)\n",
				size_col, size_row, NB_Block);

		return -1;
	}

	if (NB_Block % size_row != 0 || NB_Block % size_col != 0) {
		info_print(0, row + group_row * size_col, col + group_col * size_row,
				"The number of Block by processor is not an integer\n");
		return -1;
	}

	if (row >= size_col || col >= size_row) {
		info_print(0, row, col,
				"I'm useless bye!!! col: %zu row: %zu, " "size_col: %zu , size_row: %zu \n",
				col, row, size_col, size_row);
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

	debug_print(1, row, col,
			"I'm initialized col: %zu row: %zu, " "size_col: %zu , size_row: %zu, my rank: %d \n",
			col, row, size_col, size_row, myrank);

	/*------------------------Initialize the matrices--------------------------*/

	/* think about a common interface
	 *  int pdgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans,
	 *                 m, n, k, alpha, a, ia, ja, lda, b, ib, jb, ldb,
	 *                 beta, c, ldc, Comm, rank );
	 */

	/*------------------------Prepare the Communication Layer------------------*/
	/* Split comm size_to row and column comms */
	MPI_Comm row_comm, col_comm, group_col_comm, group_row_comm;
	MPI_Comm my_last_world;
	MPI_Comm_split(my_world, 0, key, &my_last_world);
	MPI_Comm_free(&my_world);

	MPI_Comm_split(my_last_world, group, myrank, &group_comm);

	MPI_Comm_split(my_last_world,
			col + row * size_row + size_col * size_row * group_row, group_col,
			&group_row_comm);
	MPI_Comm_split(my_last_world,
			col + row * size_row + size_col * size_row * group_col, group_row,
			&group_col_comm);
	/* color by row, rank by column */
	MPI_Comm_split(group_comm, row, col, &row_comm);
	/* color by column, rank by row */
	MPI_Comm_split(group_comm, col, row, &col_comm);
	/*-------------------------Communication Layer can be used------------------*/

	debug_print(1, row + group_row * size_col, col + group_col * size_row,
			"row %zu col %zu group_col %zu group_row %zu " "size_group_col %zu, size_group_row %zu " "group_col_comm %zu group_row_comm %zu\n",
			row, col, group_col, group_row, size_group_col, size_group_row,
			col + row * size_col + size_col * size_row * group_col,
			col + row * size_col + size_col * size_row * group_row);

	size_t ldb = n;
	size_t lda = k_a;
	MPI_Datatype Block_a;
	MPI_Datatype Block_a_local;
	MPI_Datatype Block_b;

	if (is_cyclic) {
		MPI_Type_vector(m, Block_size_group / size_row, lda, MPI_DOUBLE,
				&Block_a);

		MPI_Type_vector(m, Block_size_group / size_row, Block_size_group,
				MPI_DOUBLE, &Block_a_local);

		MPI_Type_vector(Block_size_group / size_col, n, ldb, MPI_DOUBLE,
				&Block_b);

	} else {
		MPI_Type_vector(m, Block_size_group, lda, MPI_DOUBLE, &Block_a);

		MPI_Type_vector(m, Block_size_group, Block_size_group, MPI_DOUBLE,
				&Block_a_local);

		MPI_Type_vector(Block_size_group, n, ldb, MPI_DOUBLE, &Block_b);
	}

	MPI_Type_commit(&Block_a);
	MPI_Type_commit(&Block_a_local);
	MPI_Type_commit(&Block_b);

	double * B_a;
	double * B_b;
	double * a_local;
	double * b_local;

	matrices_initialisation(&a, &b, &c, m, k_a, k_b, n,
			row + group_row * size_col, col + group_col * size_row);
	blocks_initialisation(&a_local, &b_local, m, Block_size_group, n);

	B_proc_col = k_b / Block_size_group; // Number of block on one processor
	B_proc_row = k_a / Block_size_group; // Number of block on one processor
	B_group_col = B_proc_col * size_col; // Number of block on one group
	B_group_row = B_proc_row * size_row; // Number of block on one group

	/* for each group_block broadcast it between the groups and compute it */
	MPI_Barrier(my_last_world);
	get_time(&start_time);

	/*Communicate the block between the groups */
	int iter;
	int NB_block_group = k / Block_size_group;

	/*--------------------Allocation of matrices block-------------------------*/
	double *a_Summa, *b_Summa;
	blocks_initialisation(&a_Summa, &b_Summa, m, Block_size, n);

	/*--------------------Communication types for MPI--------------------------*/
	MPI_Datatype *Block_Summa_a;
	MPI_Datatype Block_Summa_a_group_local;
	MPI_Datatype Block_Summa_a_group_global;
	MPI_Datatype Block_Summa_a_local;
	MPI_Datatype Block_Summa_b;

	MPI_Type_vector(m, Block_size, Block_size_group, MPI_DOUBLE,
			&Block_Summa_a_group_local);
	MPI_Type_vector(m, Block_size, k_a, MPI_DOUBLE,
			&Block_Summa_a_group_global);
	MPI_Type_vector(m, Block_size, Block_size, MPI_DOUBLE,
			&Block_Summa_a_local);
	MPI_Type_vector(Block_size, n, n, MPI_DOUBLE, &Block_Summa_b);

	MPI_Type_commit(&Block_Summa_a_group_local);
	MPI_Type_commit(&Block_Summa_a_group_global);
	MPI_Type_commit(&Block_Summa_a_local);
	MPI_Type_commit(&Block_Summa_b);
	/*-------------Communication types for MPI are configured-------------------*/

	for (iter = 0; iter < NB_block_group; iter++) {
		// pivot on group layer
		if (is_cyclic) {
//#ifdef CYCLIC
			//size_row * size_col * size_group_col * size_group_row
			debug_print(1, row + group_row * size_col,
					col + group_col * size_row,
					"iter : %d, B_group_row %zu, B_group_col %zu, " "B_proc_row %zu, B_proc_col %zu \n",
					iter, B_group_row, B_group_col, B_proc_row, B_proc_col);
			pivot_group_col = (iter / B_group_row);
			pivot_group_row = (iter / B_group_col);

			// Which data to send????
			// Which processor send?
			// All number of block between group should be proportional to the
			// number of processor on one row or col of a group
			int size_block_a = Block_size_group / size_row;
			int size_block_b = Block_size_group / size_col;

			debug_print(1, row + group_row * size_col,
					col + group_col * size_row,
					"pivot_group_row : %zu pivot_group_col : %zu " "pivot_row : %zu pivot_col : %zu\n",
					pivot_group_row, pivot_group_col, row, col);

			//position of the block
			start = iter * Block_size_group / Block_size;
			end = (iter + 1) * Block_size_group / Block_size;
			pos_a = (iter % B_group_row) * size_block_a;
			pos_b = (iter % B_group_col) * size_block_b * ldb;

			debug_print(1, row + group_row * size_col,
					col + group_col * size_row,
					"pos_a : %zu B_group_row : %zu " "B_proc_row : %zu B_proc_row : %zu Block_size_group: %zu\n",
					pos_a, B_group_row, B_proc_row, B_proc_row,
					Block_size_group);
		} else {
//#else
			debug_print(1, row + group_row * size_col,
					col + group_col * size_row,
					"iter : %d, B_group_row %zu, B_group_col %zu, " "B_proc_row %zu, B_proc_col %zu \n",
					iter, B_group_row, B_group_col, B_proc_row, B_proc_col);
			pivot_group_col = (iter / B_group_row);
			pivot_group_row = (iter / B_group_col);

			// Which data to send????
			// Which processor send?
			pivot_col = (iter % B_group_row) / B_proc_row;
			pivot_row = (iter % B_group_col) / B_proc_col;

			debug_print(1, row + group_row * size_col,
					col + group_col * size_row,
					"pivot_group_row : %zu pivot_group_col : %zu " "pivot_row : %zu pivot_col : %zu\n",
					pivot_group_row, pivot_group_col, pivot_row, pivot_col);

			//position of the block
			start = iter * Block_size_group / Block_size;
			end = (iter + 1) * Block_size_group / Block_size;
			pos_a = ((iter % B_group_row) % B_proc_row) * Block_size_group;
			pos_b = ((iter % B_group_col) % B_proc_col) * ldb * Block_size_group;

			debug_print(1, row + group_row * size_col,
					col + group_col * size_row,
					"pos_a : %zu B_group_row : %zu " "B_proc_row : %zu B_proc_row : %zu Block_size_group: %zu\n",
					pos_a, B_group_row, B_proc_row, B_proc_row,
					Block_size_group);
		}
//#endif

		get_time(&start_time_intern);
		//Broadcast the row
		size_t lda_local;

		if (size_group_row > 1) {
			MPI_Datatype * Block;
//#ifndef CYCLIC
			if ((!is_cyclic) && (pivot_col != col)) {
				B_a = NULL;
				lda_local = 0;
				//TODO see if it works without
				//useless but to avoid a segfault when we call Summa
				Block_Summa_a = &Block_Summa_a_group_local;
				debug_print(1, row + group_row * size_col,
						col + group_col * size_row,
						"NULL Bcast for me col %zu pivot_col %zu\n", col,
						pivot_col);
			} else {
//#endif
				if (pivot_group_col != group_col) {
					B_a = a_local;
					lda_local = Block_size_group;
					debug_print(1, row + group_row * size_col,
							col + group_col * size_row,
							"recieve group B_a col %zu pivot_col %zu\n",
							col, is_cyclic?col:pivot_col
							);
					Block = &Block_a_local;
					Block_Summa_a = &Block_Summa_a_group_local;
				} else {
					B_a = a + pos_a;
					lda_local = lda;
					debug_print(1, row + group_row * size_col,
							col + group_col * size_row,
							"sent group B_a col %zu pivot_col %zu\n",
							col, is_cyclic?col:pivot_col
							);

					Block = &Block_a;
					Block_Summa_a = &Block_Summa_a_group_global;
				}

				hpnla_bcast(B_a, 1, *Block, pivot_group_col, group_row_comm,
						bcast_algorithm);
#if(2 <= DEBUG)
				check_reception_a(a_local, m, B_k, row, col);
#endif
//#ifndef CYCLIC
			}
//#endif
		} else {
			B_a = a + pos_a;
			lda_local = lda;
			debug_print(1, row + group_row * size_col,
					col + group_col * size_row,
					"No group Bcast position of B_a size_group_row: %zu \n",
					size_group_row);
			Block_Summa_a = &Block_Summa_a_group_local;
		}

		//Broadcast the col
		if (size_group_col > 1) {

//#ifndef CYCLIC
			if ( (!is_cyclic) && (pivot_row != row) ) {
				B_b = NULL;
				debug_print(1, row + group_row * size_col,
						col + group_col * size_row,
						"NULL group B_b pos_b %zu\n", pos_b);
			} else {
				if (pivot_group_row != group_row) {
					B_b = b_local;
					debug_print(1, row + group_row * size_col,
							col + group_col * size_row,
							"Recieve group B_b pivot_row %zu row %zu\n",
							row, is_cyclic?row:pivot_row
							);
				} else {
					B_b = b + pos_b;
					debug_print(1, row + group_row * size_col,
							col + group_col * size_row,
							"sent group B_b pivot_row %zu row %zu\n",
							row, is_cyclic?row:pivot_row
							);
				}

				hpnla_bcast(B_b, 1, Block_b, pivot_group_row, group_col_comm,
						bcast_algorithm);

#if(2 <= DEBUG)
				check_reception_b(b_local, B_k, n, row, col);
#endif
//#ifndef CYCLIC
			}
//#endif
		} else {
			B_b = b + pos_b;
			debug_print(1, row + group_row * size_col,
					col + group_col * size_row,
					"No group Bcast position of B_b size_group_col: %zu \n",
					size_group_col);
		}

		MPI_Barrier(my_last_world);

		get_time(&end_time_intern);
		communication_time += get_timediff(&start_time_intern, &end_time_intern);

		debug_print(1, row + group_row * size_col, col + group_col * size_row,
				"lda_local %zu, ldb %zu, ldc %zu, " "m %zu,  k_a %zu,  k_b %zu,  n %zu, Block_size %zu, " "start %zu, end %zu, row %zu,  col %zu, " "size_row %zu,  size_col %zu\n",
				lda_local, ldb, n, m, k_a, k_b, n, Block_size, start, end, row,
				col, size_row, size_col);

//#ifndef CYCLIC
		if (is_cyclic || (pivot_row == row)) {
//#endif
			debug_print(1, row + group_row * size_col,
					col + group_col * size_row,
					"pos_a : %zu b : %e pos_b : %zu c : %e iter : %d\n", pos_a,
					B_b[0], pos_b, c[n + 1], iter);
//#ifndef CYCLIC
	//	if (is_cyclic || (pivot_col == col))
//#endif
			debug_print(1, row + group_row * size_col,
					col + group_col * size_row,
					"a : %e pos_a : %zu pos_b : %zu c : %e iter : %d\n", B_a[0],
					pos_a, pos_b, c[n + 1], iter);
		}


		debug_print(1, row + group_row * size_col, col + group_col * size_row,
				"a : %e pos_a : %zu pos_b : %zu c : %e iter : %d start : %zu end : %zu\n",
				B_a[0], pos_a, pos_b, c[1], iter, start, end);

		/* --------- Summa algorithm modify to take the block of matrices ------*/
		Summa(gemm, B_a, B_b, c, lda_local, ldb, n, m, k_a, k_b, n, Block_size,
				start, end, row, col, size_row, size_col, a_Summa, b_Summa,
				*Block_Summa_a, Block_Summa_a_local, Block_Summa_b, row_comm,
				col_comm, 1, &summa_computation_t, &summa_communication_t,
				bcast_algorithm,
				distribution);

		debug_print(1, row + group_row * size_col, col + group_col * size_row,
				"output: a : %e pos_a : %zu pos_b : %zu c : %e iter : %d start : %zu end : %zu\n",
				B_a[0], pos_a, pos_b, c[1], iter, start, end);
		/*-------------------------End Summa algorihtm--------------------------*/
	}

	/* the work is finished*/
	MPI_Barrier(my_last_world);

	get_time(&end_time);
	time = get_timediff(&start_time, &end_time);

	printf(
			"Proc: (%zu, %zu), communication time: %le, total time: %le, summa_computation_time: %le, summa_communication_time: %le, global_m: %zu, global_n: %zu, global_k: %zu, m: %zu,  k_a: %zu,  k_b: %zu,  n: %zu, Block_size: %zu, Block_size_group: %zu, size_row: %zu, size_col: %zu, size_group_row: %zu, size_group_col: %zu, bcast_alg: %zu\n",
			row + group_row * size_col, col + group_col * size_row,
			communication_time, time, summa_computation_t,
			summa_communication_t, global_m, global_n, global_k, m, k_a, k_b, n,
			Block_size, Block_size_group, size_row, size_col, size_group_row,
			size_group_col, bcast_algorithm);

	MPI_Barrier(my_last_world);

//#if CHECK_25D
	debug_print(1, row + group_row * size_col, col + group_col * size_row,
			"finished with c : %e iter : %d\n", c[m], iter);
	check_result(c, a, b, m, n, k_a, k_b, row, col, size_row * size_group_row,
			size_col * size_group_col);
//#endif

	// close properly the pragram
	MPI_Type_free(&Block_Summa_a_group_local);
	MPI_Type_free(&Block_Summa_a_group_global);
	MPI_Type_free(&Block_Summa_a_local);
	MPI_Type_free(&Block_Summa_b);
	MPI_Type_free(&Block_b);
	MPI_Type_free(&Block_a_local);
	MPI_Type_free(&Block_a);

	SMPI_SHARED_FREE(a_Summa);
	SMPI_SHARED_FREE(b_Summa);

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

	return 0;
}
