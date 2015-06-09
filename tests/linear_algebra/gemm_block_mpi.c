/*!
 * Block Matrix Multiplication example
 *
 * Authors: Quintin Jean-Noël
 */

#include "config.h"
#include "tools/hpnla_debug.h"
#include "Matrix_init.h"
#include "Summa.h"
#include "hpnla_cblas.h"
#include "tools/hpnla_timer.h"

#include <mpi.h>
#include <math.h>
#include <getopt.h>
#include <stdio.h>
#include <stdbool.h>
#include <sys/time.h>
#ifdef HDNLA_SMPI
#include <smpi.h>
#else
#define SMPI_SHARED_MALLOC malloc
#define SMPI_SHARED_FREE free
#endif

int main(int argc, char ** argv) {
	double *a, *b, *c; //matrices
	double summa_computation_t = 0, summa_communication_t = 0;

	size_t m = 1024, n = 1024;
	size_t k_a = 1024, k_b = 1024; // matrix a, b local k sizes
	size_t NB_Block = 16, Block_size = k_a / NB_Block;
	/* x index on M
	 y index on N
	 Z index on K */

	int bcast_algorithm = 5; //original algorithm inside mpi
	int distribution = 0;    // default matrix distribution is non-cyclic
	int myrank;
	int NB_proc;
	size_t row, col, size_row, size_col; //description: vitual processor topology
	row = 0;
	col = 0;

	MPI_Init(&argc, &argv);

	/* Find out my identity in the default communicator */

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &NB_proc);

	if (NB_proc != 1)
		for (size_col = NB_proc / 2; NB_proc % size_col; size_col--)
			;
	else
		size_col = 1;

	size_row = NB_proc / size_col;

	if (size_row > size_col) {
		size_col = size_row;
		size_row = NB_proc / size_col;
	}

	// for the degub
#if DEBUG_MPI
	size_t loop=1;
	while(loop==1);
#endif

	int opt;
	optind = 1;

	//get the parameter from command line
	while ((opt = getopt(argc, argv, "hr:d:a:M:N:K:B:")) != -1) {
		switch (opt) {
		case 'h':
			info_print(0, row, col,
					"Usage: mxm_cblas_test [options]\n" ""
							"	-M I	M size (default: %zu)\n"
							"	-N I	N size (default: %zu)\n"
							"	-K I	K size (default: %zu)\n"
					        "   -B I	Block size on the k dimension(default: %zu)\n"
					        "   -r I	processor row size (default: %zu)\n"
					        "   -c I	processor col size (default: %zu)\n"
					        "   -a I    bcast algorithm (default: %d)\n"
					        "   -d I    is_cyclic distribution or not (default: %d)\n"
					        "   -h	help\n",
					m, n, k_a, Block_size, size_row, size_col, bcast_algorithm,
					distribution);
			return 0;
		case 'M':
			m = atoi(optarg);
			break;
		case 'N':
			n = atoi(optarg);
			break;
		case 'K':
			k_a = atoi(optarg);
			k_b = k_a;
			break;
		case 'B':
			Block_size = atoi(optarg);
			break;
		case 'r':
			size_row = atoi(optarg);
			break;
		case 'c':
			size_col = atoi(optarg);
			break;
		case 'a':
			bcast_algorithm = atoi(optarg);
			break;
		case 'd':
			distribution = atoi(optarg);
			break;
		}
	}

	if (NB_proc < (int) (size_row * size_col)) { //No of process should be equal to No. of process in row and column
		info_print(0, row, col,
				"Not enough processors NB_proc : %d required : %zu\n", NB_proc,
				size_row * size_col);
		exit(-1);
	}

	NB_proc = size_row * size_col; // Check in case when user enter more processes than required
								   // then use only required & release other unnecessary processes
	row = myrank / size_row;
	col = myrank % size_row;

	/*-------------------------Check some mandatory conditions-----------------*/
	NB_Block = k_a / Block_size;
	if (k_a % Block_size != 0) {
		info_print(0, row, col,
				"The matrix size has to be proportionnal to" " the number of blocks: %zu\n",
				NB_Block);
		exit(-1);
	}

	if (size_row > NB_Block || size_col > NB_Block) {
		info_print(0, row, col,
				"Number of blocks is too small compare to " "the number of processors (%zu,%zu) in a row or a col (%zu)\n",
				size_col, size_row, NB_Block);
		exit(-1);
	}

	if (NB_Block % size_row != 0 || NB_Block % size_col != 0) {
		info_print(0, row, col, "The number of Block by processor is not an %s",
				"size_terger\n");
		exit(-1);
	}
	if (row >= size_col || col >= size_row) {
		info_print(0, row, col,
				"I'm useless bye!!! col: %zu row: %zu, " "size_col: %zu , size_row: %zu \n",
				col, row, size_col, size_row);

		/*----------------------Prepare the Communication Layer------------------*/
		/* add useless processor on a new color to execute the matrix
		 * multiplication with the other processors*/

		/* Split comm size_to row and column comms */
		MPI_Comm row_comm, col_comm;
		/* color by row, rank by column */
		MPI_Comm_split(MPI_COMM_WORLD, size_row, col, &row_comm);
		/* color by column, rank by row */
		MPI_Comm_split(MPI_COMM_WORLD, size_col, row, &col_comm);
		/*------------------------Communication Layer can be used----------------*/

		MPI_Barrier(MPI_COMM_WORLD );
		exit(0);
	}
	info_print(0, row, col,
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
	MPI_Comm row_comm, col_comm;
	/* color by row, rank by column */
	MPI_Comm_split(MPI_COMM_WORLD, row, col, &row_comm);
	/* color by column, rank by row */
	MPI_Comm_split(MPI_COMM_WORLD, col, row, &col_comm);
	/*-------------------------Communication Layer can be used-----------------*/

	// matrix sizes
	m = m / size_col;
	n = n / size_row;
	k_a = k_a / size_row;
	k_b = k_b / size_col;

	// Defined the device if we use the GPU
	//TODO explain parameters
	hpnla_gemm * gemm = hpnla_gemm_alloc(NULL );

	matrices_initialisation(&a, &b, &c, m, k_a, k_b, n, row, col);

	debug_print(1, row, col,
			"m %zu,  k_a %zu,  k_b %zu,  n %zu,  Block_size %zu, " "group*NB_Block/NB_groups %zu, (group+1)*NB_Block/NB_groups %zu," " row %zu,  col %zu,  size_row %zu,  size_col %zu\n",
			m, k_a, k_b, n, Block_size, (size_t )0, NB_Block, row, col,
			size_row, size_col);
	double time;
	struct timespec start_time, end_time; //time mesure
	get_time(&start_time);
	//  Summa(gemm, a, b, c, m, k_a, k_b, n, Block_size, 0, NB_Block,
	//         row, col, size_row, size_col, row_comm, col_comm);
	/*-------------------Configuration for Summa algorihtm--------------------*/
	/*--------------------Allocation of matrices block-------------------------*/
	double *a_Summa, *b_Summa;
	blocks_initialisation(&a_Summa, &b_Summa, m, Block_size, n);

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

	Summa(gemm, a, b, c, k_a, n, n, m, k_a, k_b, n, Block_size, 0, NB_Block,
			row, col, size_row, size_col, a_Summa, b_Summa, Block_a,
			Block_a_local, Block_b, row_comm, col_comm, 0, &summa_computation_t,
			&summa_communication_t, bcast_algorithm, distribution);

	MPI_Type_free(&Block_a);
	MPI_Type_free(&Block_a_local);
	MPI_Type_free(&Block_b);

	SMPI_SHARED_FREE(a_Summa);
	SMPI_SHARED_FREE(b_Summa);
	/*-------------------------End Summa algorihtm----------------------------*/
	get_time(&end_time);
	time = get_timediff(&start_time, &end_time);
	printf("time: %le nanoseconds\n", time);

	if (size_row * size_col == 1) {

		/* inital code for dgemm */
		size_t lda = k_a, ldb = n, ldc = n; // matrix line size
		double * check = (double *) SMPI_SHARED_MALLOC(sizeof(double) * m * n);
		if (check == 0) {
			perror("Error allocation Matrix Check");
			exit(-1);
		}
		size_t x, y;
		for (x = 0; x < m; x++) {
			for (y = 0; y < n; y++) {
				check[x * ldc + y] = (double) 0;
			}
		}

		hpnla_gemm_execute(
				gemm, //user parameters: here it's equal to NULL
				CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k_a, 1, a, lda,
				b, ldb, 1, check, ldc);

		/*Display for checking */bool result_good = true;
		for (x = 0; x < m; x++) {
			for (y = 0; y < n; y++) {
				/* WARNING this could be lead to some errors (precision with double)*/
				if (c[x * ldc + y] != check[x * ldc + y]) {
					result_good = false;
					printf("%lf\t%lf\n", c[x * ldc + y], check[x * ldc + y]);
				}
			}
		}
		NOT_USED(result_good);
	}

	check_result(c, a, b, m, n, k_a, k_b, row, col, size_row, size_col);

	SMPI_SHARED_FREE(a);
	SMPI_SHARED_FREE(b);
	SMPI_SHARED_FREE(c);

	hpnla_gemm_free(gemm);

	MPI_Barrier(MPI_COMM_WORLD );
	MPI_Comm_free(&row_comm);
	MPI_Comm_free(&col_comm);

        MPI_Finalize();

	return 0;
}
