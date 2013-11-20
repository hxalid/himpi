/*
 * Block Matrix Multiplication example
 *
 * Authors: Quintin Jean-Noël
 */


#include "config.h"
#include "tools/hpnla_debug.h"
#include "tools/hpnla_bind.h"
#include "tools/hpnla_conf.h"
#include "Matrix_init.h"
#include "2.5D_MM.h"
#include "cblas_wrappers/hpnla_cblas.h"

/*int sched_setaffinity(pid_t pid, size_t cpusetsize, cpu_set_t *mask);
  int sched_getaffinity(pid_t pid, size_t cpusetsize, cpu_set_t *mask);
 */
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

  size_t    m   = 1024 , n   = 1024 , k = 1024;
  size_t    NB_Block = 16;
  size_t    Block_size = k/NB_Block ;
  size_t    NB_groups = 1, group = 0, key = 0;
  /* x index on M
     y index on N
     Z index on K */


  MPI_Comm comm = MPI_COMM_WORLD;
  int root = 0;
  int bcast_algorithm = 4; //original alg inside mpi
  int distribution = 0;    // default matrix distribution is non-cyclic
  int myrank;
  int NB_proc;
  size_t row, col, size_row, size_col; //description: vitual processor topology
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
  char *conf_file = "machinefile";
  optind = 1;

  //get the parameter from command line
  while ((opt = getopt(argc, argv, "hdf:r:c:s:a:d:M:N:K:B:G:g:k:P:")) != -1) {
    switch(opt) {
      case 'h':
        info_print(0, row, col,
            "Usage: mxm_cblas_test [options]\n"
            "	-M I	M size (default: %zu)\n"
            "	-N I	N size (default: %zu)\n"
            "	-K I	K size (default: %zu)\n"
            "	-B I	Block size on the k dimension(default: %zu)\n"
            "	-G I	Number of processor groups(default: %zu)\n"
            "	-g I	group index(default: %zu)\n"
            "	-k I	group rank(default: %zu)\n"
            "	-r I	processor row size (default: %zu)\n"
            "	-c I	processor col size (default: %zu)\n"
            "   -a I    broadcast algorithm (default: %zu)\n"
            "   -f {Filename} provide the file with the configuration\n"
            "	-s  display the configuration file on the stderr\n"
        	"   -d I    is_cyclic distribution or not (default: %d)\n"
            "	-h	help\n",
            m, n, k, Block_size, NB_groups, group, key, row, col, bcast_algorithm);
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
      case 'G':
        NB_groups = atoi(optarg);
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
    hpnla_print_conf(MPI_COMM_WORLD, root, stdout, "cpu", "");
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(0);
  }

  hpnla_process_conf conf = hpnla_get_conf(MPI_COMM_WORLD, conf_file);

  if (hpnla_bind_process(conf.bind) == -1) goto end;


    //    info_print(0, row, col, "group the process: %s\n",conf[4]);
    if(NB_groups !=1){
      /* Get my group number from the config file */
      group = (size_t)atoi(conf.subopts);
    }

    /*
    if(conf.argv[2] != NULL){
      myrank = (size_t)atoi(conf.argv[2]);
      info_print(0, row, col, "process the renumarotation\n");
      MPI_Comm_split(MPI_COMM_WORLD, 0, myrank, &comm);
      int tmp_rank;
      MPI_Comm_rank( comm, &tmp_rank );
      if(tmp_rank != myrank)
        info_print(0, row, col, "error during the renumerotation got %d instead of %d\n",tmp_rank,myrank);
    }*/
  
  // Defined the device if we use the GPU
  //TODO explain parameters
  hpnla_gemm  * gemm  = hpnla_gemm_alloc(&conf);

  two_dot_five(gemm, m, k, n, Block_size, group, key,
      size_row, size_col,  NB_groups, comm, bcast_algorithm, distribution);

  // close properly the pragram
  hpnla_gemm_free(gemm);
end:
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
