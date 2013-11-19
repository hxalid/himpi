#include "Matrix_init.h"
#include "tools/hdnla_debug.h"
#include "tools/hdnla_timer.h"
#include "Cannon.h"
#include <sys/time.h>
#include <stdbool.h>


double cannon(hdnla_gemm  * gemm,
		double *a, double *b, double *c,
		size_t local_matrix_dim,int proc_assign, int *cart_coords, MPI_Comm comm_cart)
{
	size_t err;
	int up_rank, down_rank, left_rank, right_rank;
	int shift_src, shift_dest;
	double alpha = 1, beta = 1;  //C := alpha * a * b + beta * c
	size_t lda = local_matrix_dim, ldb = local_matrix_dim,ldc = local_matrix_dim;
	MPI_Status status;

	init_timer();
	double communication_time = 0, shift_time = 0;
	struct timespec start_time_intern, end_time_intern; //time mesure for Cannon without considering initial alligment time
	struct timespec shift_stime, shift_etime; // time measuse for shift inside cannon computation


	/*************************Distributed Matrix Multiplication algorithm**********************************/

	MPI_Cart_shift(comm_cart, 1, -1, &right_rank, &left_rank); //Compute rank of left shift processor
	MPI_Cart_shift(comm_cart, 0, -1, &down_rank, &up_rank);   // Compute rank of up shift processor 

	//Perform initial matrix alignment for A 
	MPI_Cart_shift(comm_cart, 1, -cart_coords[0], &shift_src, &shift_dest); 
	//printf("Rank: %d, Coordinates:(%d, %d),source rank: %d, Destination rank:%d\n"
	//       ,cart_rank, cart_coords[0], cart_coords[1],shift_src,shift_dest);   
	MPI_Sendrecv_replace(a, local_matrix_dim *local_matrix_dim,  MPI_DOUBLE, shift_dest,
			1, shift_src, 1, comm_cart, &status); 

	//Perform initial matrix alignment for B 
	MPI_Cart_shift(comm_cart, 0, -cart_coords[1], &shift_src, &shift_dest); 
	MPI_Sendrecv_replace(b, local_matrix_dim *local_matrix_dim , MPI_DOUBLE, shift_dest, 
			1, shift_src, 1, comm_cart, &status); 

	err = MPI_Barrier(comm_cart);

	get_time(&start_time_intern);

	size_t iter;
	for( iter = 0; iter < (size_t)proc_assign; iter++ ){ // start =0 end =No_of_process in row-1

		hdnla_gemm_execute(gemm, //user parameters: here it's equal to NULL
				CblasRowMajor, CblasNoTrans, CblasNoTrans,local_matrix_dim, 
				local_matrix_dim,local_matrix_dim,alpha, a, lda, b, ldb, beta, c, ldc);


      /********************************* Shift matrix a left by one ******************************************/ 
		get_time(&shift_stime);
		if (iter == (size_t)proc_assign-1) goto end;
		err =MPI_Sendrecv_replace(a, local_matrix_dim *local_matrix_dim , MPI_DOUBLE, 
				left_rank, 1, right_rank, 1, comm_cart, &status); 

		if (err != MPI_SUCCESS) {
			perror("Error in Shift B\n");
			exit(-1);
		}

      /************************************ Shift matrix b up by one ******************************************/ 
		err=MPI_Sendrecv_replace(b, local_matrix_dim *local_matrix_dim , MPI_DOUBLE, 
				up_rank, 1, down_rank, 1, comm_cart, &status); 

		if (err != MPI_SUCCESS) {
			perror("Error in Shift B\n");
			exit(-1);

		}
           end:
        	get_time(&shift_etime);
		shift_time += get_timediff( &shift_stime, &shift_etime);

	}

  info_print(0,(size_t)cart_coords[0],(size_t)cart_coords[1], "Cannon inside shift time: %le microseconds\n",shift_time);
	get_time(&end_time_intern);

	communication_time = get_timediff(&start_time_intern,&end_time_intern);
  info_print(0, (size_t)cart_coords[0], (size_t)cart_coords[1], "Cannon inside computation_time: %le microseconds\n"
		    ,communication_time);


	/****************************************Check Result *************************************************/
	check_result(c, a, b, local_matrix_dim, local_matrix_dim,local_matrix_dim,local_matrix_dim,
		      cart_coords[0], cart_coords[1],proc_assign, proc_assign);
	return (0);
}
