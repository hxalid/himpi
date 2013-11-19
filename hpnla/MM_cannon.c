#include "config.h"
#include "tools/hdnla_debug.h"
#include "Matrix_init.h"
#include "hdnla_cblas.h"
#include "hdnla_memory.h"
#include "tools/hdnla_timer.h"
# include "Cannon.h"

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


int main(int argc, char ** argv)
{
  double *a  , *b  , *c  ;     // Matrices on each processor
  size_t matrix_dim=256;      // Matrices dimensions
  size_t local_matrix_dim;    // Determine the dimension of the Local matrices blocks
  int myrank;
  int NB_proc;
  int proc_assign;

  // for the degub
  #if DEBUG_MPI
  size_t loop=1;
  while(loop==1);
   #endif

  int opt;
  optind = 1;

  //get the parameter from command line
  while ((opt = getopt(argc, argv, "hM:")) != -1) {
    switch(opt) {
  
      case 'h':
        info_print(0, (size_t)0,(size_t) 0,
                    "Usage: mxm_cblas_test [options]\n"
                    "	-M I	Matrix Dimensions (default: %zu)\n"
                    "	-h	help\n",
                    (size_t)matrix_dim);
        return 0;
      case 'M':
        matrix_dim = atoi(optarg);
        break;
    }
  }

  int err;
	NOT_USED(err);
  err= MPI_Init(&argc, &argv);

  /* Find out my identity in the default communicator */

  err=MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  err=MPI_Comm_size (MPI_COMM_WORLD, &NB_proc);

  proc_assign      = sqrt(NB_proc);
  local_matrix_dim = matrix_dim/ proc_assign;


 /************************************ Set up the Cartesian Topology for processors***************************************/
  int dims [2], periods [2];
  int cart_rank;
  int cart_coords[2]; 
  MPI_Comm comm_cart; 

  dims [0] = dims [1]= proc_assign; //2-D Topology
  periods [0]= periods [1]= 1; //It indicate whether the dimension is cyclic or not, 1 mean both dimensions are cyclic
  MPI_Cart_create (MPI_COMM_WORLD,2,dims,periods,1, &comm_cart); //Create new communicator for 2-D cartesion topology with rank reordering
  MPI_Comm_rank(comm_cart, & cart_rank); //Get rank of processors in new Topology
  MPI_Cart_coords(comm_cart, cart_rank, 2,  cart_coords);// Get the Coordinates of processors in new Topology
 
  info_print(0,(size_t) cart_coords[0],(size_t) cart_coords[1], "I'm initialized Row: %d Col: %d, "
              "Size_row: %d , Size_col: %d, My rank: %d, Matrix Dimention: %zu \n"
              ,cart_coords[0],cart_coords[1],proc_assign,proc_assign, cart_rank, matrix_dim);


  // Defined the device if we use the GPU
  //TODO explain parameters
  hdnla_gemm  * gemm  = hdnla_gemm_alloc(NULL);

/********************************************************************************************************************/

  matrices_initialisation (&a, &b, &c, local_matrix_dim, local_matrix_dim ,local_matrix_dim,
		            local_matrix_dim,cart_coords[0], cart_coords[1]); 

  double time;
  struct timespec start_time, end_time; //time mesure

  get_time(&start_time);
  cannon(gemm, a, b, c, local_matrix_dim, proc_assign,cart_coords, comm_cart);
  get_time(&end_time);

  time = get_timediff(&start_time,&end_time);
  info_print(0, (size_t)cart_coords[0], (size_t)cart_coords[1], "Cannon computation + communication Time: %le microseconds\n",time);

  err = MPI_Barrier(comm_cart);
 
  SMPI_SHARED_FREE(a);
  SMPI_SHARED_FREE(b);
  SMPI_SHARED_FREE(c);

  MPI_Comm_free(&comm_cart); /* Free the communicator */
  hdnla_gemm_free(gemm);

  MPI_Finalize(); 

  return 0;
}


