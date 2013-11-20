/*!
 * Classical Block Matrix Multiplication example
 *
 * Authors: Quintin Jean-Noël
 */
#include "config.h"
#include "Matrix_init.h"
#include "tools/hpnla_debug.h"
#include "tools/hpnla_timer.h"
#include "communication/hpnla_bcast.h"
#include "cblas_wrappers/hpnla_cblas.h"
#include "Summa.h"
#include <sys/time.h>
#ifdef HDNLA_SMPI
#include <smpi.h>
#else
#define SMPI_SAMPLE_GLOBAL(x,y) do{}while(0);
#endif

inline double Summa(hpnla_gemm * gemm, 
		    double *a, double *b, double *c,
		size_t lda, size_t ldb, size_t ldc, size_t m, size_t k_a, size_t k_b,
		size_t n, size_t Block_size, size_t start, size_t end, size_t row,
		size_t col, size_t size_row, size_t size_col, double *a_local,
		double *b_local, MPI_Datatype Block_a, MPI_Datatype Block_a_local,
		MPI_Datatype Block_b, MPI_Comm row_comm, MPI_Comm col_comm, int subs,
		double *computation_time, double *communication_time,
		int bcast_algorithm,
		int distribution) {

	double *B_a, *B_b; 				// matrix blocks
	double alpha = 1, beta = 1;  	// C := alpha * a * b + beta * c
	size_t B_proc_col, B_proc_row;  // Number of bloc(row or col) on one processor
	B_proc_col = k_b / Block_size;  // Number of block on one processor
	B_proc_row = k_a / Block_size;  // Number of block on one processor

	int is_cyclic = (distribution > 0);

	size_t lda_local = lda;
	size_t ldb_local = ldb;

	init_timer();

	double time;
	struct timespec start_time, end_time; 				//time mesure
	struct timespec start_time_intern, end_time_intern; //time mesure

	get_time(&start_time);

	/*-------------Distributed Matrix Multiplication algorithm-----------------*/
	size_t iter;
	for (iter = start; iter < end; iter++) {
		size_t pivot_row, pivot_col, pos_a, pos_b;

	if (is_cyclic) {
		// pivot row on processor layer
		pivot_row = (iter % size_col);
		pivot_col = (iter % size_row);

		//position of the block
		if (subs == 1) {
			pos_a = (size_t) ((iter - start) / size_row) * Block_size;
			pos_b = (size_t) ((iter - start) / size_col) * ldb * Block_size;
		} else {
			pos_a = (size_t) (iter / size_row) * Block_size;
			pos_b = (size_t) (iter / size_col) * ldb * Block_size;
		}

	} else {
		// pivot row on processor layer
		pivot_row = (size_t)(iter / B_proc_col) % size_col;
		pivot_col = (size_t)(iter / B_proc_row) % size_row;

		//position of the block
		if(subs == 1) {
			pos_a = ((iter - start) % B_proc_row) * Block_size;
			pos_b = ((iter - start) % B_proc_col) * ldb * Block_size;
		} else {
			pos_a = (iter % B_proc_row) * Block_size;
			pos_b = (iter % B_proc_col) * ldb * Block_size;
		}
	}

	debug_print(3, row, col,
				"pivot: %zu, iter: %zu, B_proc_col: %zu, " "size_col:%zu, size_row: %zu\n",
				pivot_row, iter, B_proc_row, size_col, size_row);
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

			hpnla_bcast(B_a, 1, *Block, pivot_col, row_comm, bcast_algorithm);

#if(2 <= DEBUG)
			check_reception_a(a_local, m, Block_size, row, col);
#endif
		} else {
			B_a = a + pos_a;
			debug_print(1, row, col, "position of B_a: %zu \n", pos_a);
		}

		//Broadcast the col
		if (size_col > 1) {
			if (pivot_row == row) {
				B_b = b + pos_b;
				debug_print(1, row, col, "sent B_b Block_size: %zu, pos:%zu \n",ldb, pos_b);
			} else {
				B_b = b_local;
				debug_print(1, row, col, "recieve B_b %zu,%zu \n", Block_size, n);
			}

			hpnla_bcast(B_b, 1, Block_b, pivot_row, col_comm, bcast_algorithm);


#if(2 <= DEBUG)
			check_reception_b(b_local, Block_size, n, row, col);
#endif
		} else {
			B_b = b + pos_b;
			debug_print(1, row, col, "position of B_b: %zu \n", pos_b);
		}

		get_time(&end_time_intern);
		*communication_time += get_timediff(&start_time_intern, &end_time_intern);

#if DEBUG
		MPI_Barrier(row_comm);
		MPI_Barrier(col_comm);
#endif
		get_time(&start_time_intern);

		debug_print(1, row, col, "execute Gemm number: %zu\n", iter);
		//We have recieved a line of block and a colomn
		{
			SMPI_SAMPLE_GLOBAL(30, 0.1)
			{
				hpnla_gemm_execute(
						gemm, //user parameters: here equal NULL
						CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n,
						Block_size, alpha, B_a, lda_local, B_b, ldb_local, beta,
						c, ldc);
			}
		}

		get_time(&end_time_intern);
		*computation_time += get_timediff(&start_time_intern, &end_time_intern);

	}

#if DEBUG
	MPI_Barrier(row_comm);
	MPI_Barrier(col_comm);
#endif

	get_time(&end_time);
	time = get_timediff(&start_time, &end_time);

//  printf("communication time: %le nanoseconds, "  "computation time: %le nanoseconds\n",  communication_time, computation_time);

#if 0
	CHECK_Summa
	check_result(c, a, b, m, n, k_a, k_b, row, col, size_row, size_col);
#endif

	return time;
}
