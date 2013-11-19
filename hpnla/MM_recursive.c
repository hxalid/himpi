/*
 * Block Matrix Multiplication example
 *
 * Authors: Quintin Jean-Noël
 */


#include "config.h"
#include "tools/hpnla_conf.h"
#include "tools/hpnla_debug.h"
#include "tools/hpnla_bind.h"
#include "Matrix_init.h"
#include "Recursive.h"
#include "hpnla_cblas.h"
#include "hpnla_memory.h"

#include <mpi.h>
#include <math.h>
#include <getopt.h>
#include <stdio.h>
#include <string.h>



/*the groups could be created from the call
  hpnla_comm_intra hpnla_kernel.h.
  Since each proc has a local identificator the configuration could be done.
  Why have we to be directly depending on fupermod?*/

int main(int argc, char ** argv)
{

  int bcast_algorithm = 4; //original algorithm inside mpi
  int distribution = 0;    //default matrix distribution is non-cyclic
  size_t    m   = 1024 , n   = 1024 , k = 1024;
  size_t    NB_Block = 16;
  size_t    Block_size = k/NB_Block ;
  size_t    Block_group_size = Block_size;
  size_t    size_group_col = 1, size_group_row = 1, group = 0, key = 0;


  int myrank;
  int NB_proc;
  size_t row, col, size_row, size_col;  //description: vitual processor topology
  row = 0;
  col = 0;

  MPI_Init(&argc, &argv);

  /* Find out my identity in the default communicator */

  MPI_Comm_rank ( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size ( MPI_COMM_WORLD, &NB_proc );

  if(NB_proc != 1)
    for (size_col=NB_proc/2; NB_proc%size_col; size_col--);
  else
    size_col = 1;

  size_row = NB_proc/size_col;
  if (size_row > size_col){
    size_col = size_row; 
    size_row = NB_proc/size_col;
  }


  // for the degub
#if DEBUG_MPI
  size_t loop=1;
  while(loop==1);
#endif

  int opt, display = 0;
  char *conf_file = NULL;
  optind = 1;

  //get the parameter from command line
  while ((opt = getopt(argc, argv, "hdf:r:c:a:s:d:M:N:K:T:B:R:C:g:k:P:")) != -1) {
    switch(opt) {
      case 'h':
        info_print(0, row, col,
                    "Usage: mxm_cblas_test [options]\n"
                    "	-M I	M size (default: %zu)\n"
                    "	-N I	N size (default: %zu)\n"
                    "	-K I	K size (default: %zu)\n"
                    "	-B I	Block size on the k dimension(default: %zu)\n"
                    "	-T I	Block group size on the k dimension(default: %zu)\n"
                    "	-R I	processor row size (default: %zu)\n"
                    "	-C I	processor col size (default: %zu)\n"
                    "	-g I	group index(default: %zu)\n"
                    "	-k I	group rank(default: %zu)\n"
                    "	-r I	processor row size (default: %zu)\n"
                    "	-c I	processor col size (default: %zu)\n"
		            "   -a I    broadcast algorithm (default: %d)\n"
                    "   -f {Filename} provide the file with the configuration\n"
                    "	-s  display the configuration file on the stderr\n"
        		    "   -d I    is_cyclic distribution or not (default: %d)\n"
                    "	-h	help\n",
                    m, n, k, Block_size, Block_group_size,
                    size_group_col, size_group_row, group, key, row, col, bcast_algorithm, distribution);
        return 0;
      case 'M':
        m = atoi(optarg);
        break;
      case 'N':
        n   = atoi(optarg);
        break;
      case 'K':
        k  = atoi(optarg);
        break;
      case 'B':
        Block_size = atoi(optarg);
        break;
      case 'T':
        Block_group_size = atoi(optarg);
        break;
      case 'R':
        size_group_row = atoi(optarg);
        break;
      case 'C':
        size_group_col = atoi(optarg);
        break;
      case 'g':
        group = atoi(optarg);
        break;
      case 'k':
        key = atoi(optarg);
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
      case 'f':
        conf_file = strdup(optarg);
        break;
      case 's':
        display = 1;
        break;
      case 'd':
        distribution = atoi(optarg);
        break;
    }
  }

  if( display == 1 ){
    hpnla_print_conf(MPI_COMM_WORLD, myrank, stderr, "", "");
    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);
  }

  MPI_Comm world = MPI_COMM_WORLD;
  hpnla_process_conf conf;
  if(conf_file == NULL){
    info_print(0, row, col,
                "No configuration file\n");
  }else{
    conf = hpnla_get_conf(MPI_COMM_WORLD, conf_file);
  }


 // if(conf.argv == NULL){
  if (conf.subopts[0] == '\0') {
    info_print(0, row, col,
                "No configuration for me inside the file\n");
  }else{

     info_print(0, row, col, "my_rank:%d, my_group:%s\n", conf.rank_intra, conf.subopts);

     if(hpnla_bind_process(conf.bind) == -1) goto end;


     group = (size_t)atoi(conf.subopts);

/*    if(conf.argv[4] != NULL){
      myrank = (size_t)atoi(conf.argv[4]);
      info_print(0, row, col, "process the renumarotation\n");
      MPI_Comm_split(MPI_COMM_WORLD, 0, myrank, &world);
      int tmp_rank;
      MPI_Comm_rank( world, &tmp_rank );
      if(tmp_rank != myrank)
        info_print(0, row, col, "error during the renumerotation got %d instead of %d\n",tmp_rank,myrank);
    } */

  }




  // Defined the device if we use the GPU
  //TODO explain parameters
  hpnla_gemm  * gemm  = hpnla_gemm_alloc(&conf);

  recursive(gemm,m, k, n, Block_size, Block_group_size, group, key,
      size_row, size_col, size_group_row , size_group_col, world, bcast_algorithm, distribution);

  // close properly the pragram
  hpnla_gemm_free(gemm);
  
end:
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
