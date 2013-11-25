/*
 * 2D LU factorization
 *
 * Authors: Khalid Hasanov
 */


#include "config.h"
#include "tools/hpnla_conf.h"
#include "tools/hpnla_debug.h"
#include "tools/hpnla_bind.h"
#include "Matrix_init.h"
#include "lu_hfactorization.h"
#include "hpnla_cblas.h"

#ifdef HDNLA_SMPI
#include <smpi.h>
#else
#define SMPI_SHARED_FREE free
#endif

#include <mpi.h>
#include <math.h>
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char ** argv) {

    HLU_data* hlu_data = malloc(sizeof (HLU_data));
    if (hlu_data==NULL) {
         printf("malloc failed\n");
         return -1;
    }
    Platform_data* platform_data = malloc(sizeof (Platform_data));

    hlu_data->Block_size_in = 32;
    hlu_data->Block_size_out = hlu_data->Block_size_in;
    hlu_data->m_global = 1024;
    hlu_data->n_global = 1024;
    hlu_data->k_global = 1024;
    hlu_data->bcast_algorithm = 4; //original algorithm inside mpi


    platform_data->size_group_col = 1;
    platform_data->size_group_row = 1;


   // int myrank;
  //  int NB_proc;
    size_t row, col;
    row = 0;
    col = 0;

    int opt, display = 0;
    char *conf_file = NULL;
    hpnla_process_conf conf;



    MPI_Init(&argc, &argv);

    platform_data->comm = MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &platform_data->my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &platform_data->nb_proc);

    if (platform_data->nb_proc != 1)
        for (platform_data->size_col = platform_data->nb_proc / 2; 
                platform_data->nb_proc % platform_data->size_col; platform_data->size_col--);
    else
        platform_data->size_col = 1;

    platform_data->size_row = platform_data->nb_proc / platform_data->size_col;
    if (platform_data->size_row > platform_data->size_col) {
        platform_data->size_col = platform_data->size_row;
        platform_data->size_row = platform_data->nb_proc / platform_data->size_col;
    }


    // for debug
#if DEBUG_MPI
    size_t loop = 1;
    while (loop == 1);
#endif


    optind = 1;

    //get the parameter from command line
    while ((opt = getopt(argc, argv, "hdf:r:c:a:s:d:M:N:K:T:B:R:C:g:k:P:")) != -1) {
        switch (opt) {
            case 'h':
                info_print(0, row, col,
                        "Usage: LU_factorization [options]\n"
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
                        hlu_data->m_global,
                        hlu_data->n_global,
                        hlu_data->k_global,
                        hlu_data->Block_size_in,
                        hlu_data->Block_size_out,
                        platform_data->size_group_col,
                        platform_data->size_group_row,
                        hlu_data->group,
                        hlu_data->key,
                        row,
                        col,
                        hlu_data->bcast_algorithm,
                        hlu_data->distribution);
                return 0;
            case 'M':
                hlu_data->m_global = atoi(optarg);
                break;
            case 'N':
                hlu_data->n_global = atoi(optarg);
                break;
            case 'K':
                hlu_data->k_global = atoi(optarg);
                break;
            case 'B':
                hlu_data->Block_size_in = atoi(optarg);
                break;
            case 'T':
                hlu_data->Block_size_out = atoi(optarg);
                break;
            case 'R':
                platform_data->size_group_row = atoi(optarg);
                break;
            case 'C':
                platform_data->size_group_col = atoi(optarg);
                break;
            case 'g':
                hlu_data->group = atoi(optarg);
                break;
            case 'k':
                hlu_data->key = atoi(optarg);
                break;
            case 'r':
                platform_data->size_row = atoi(optarg);
                break;
            case 'c':
                platform_data->size_col = atoi(optarg);
                break;
            case 'a':
                hlu_data->bcast_algorithm = atoi(optarg);
                break;
            case 'd':
                hlu_data->distribution = atoi(optarg);
                break;
            case 'f':
                conf_file = strdup(optarg);
                break;
            case 's':
                display = 1;
                break;
        }
    }
    
    
    platform_data->nb_requested_proc = platform_data->size_row*platform_data->size_col*platform_data->size_group_row*platform_data->size_group_col;
    

    if (display == 1) {
        hpnla_print_conf(MPI_COMM_WORLD, platform_data->my_rank, stderr, "", "");
        MPI_Barrier(MPI_COMM_WORLD);
        exit(0);
    }

    if (conf_file == NULL) {
        info_print(0, row, col,
                "No configuration file\n");
    } else {
        conf = hpnla_get_conf(MPI_COMM_WORLD, conf_file);
    }


    if (conf.subopts[0] == '\0') {
        info_print(0, row, col,
                "No configuration for me inside the file\n");
    } else {

      //  info_print(0, row, col, "my_rank:%d, my_group:%s\n", conf.rank_intra, conf.subopts);

        if (hpnla_bind_process(conf.bind) == -1) goto end;

        hlu_data->group = (size_t) atoi(conf.subopts);
    }

    hlu_data->gemm = hpnla_gemm_alloc(&conf);

    
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
  /*  printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i)
        sleep(5);
    */
    
    lu_hfactorize(hlu_data, platform_data);

    // close the resources
    hpnla_gemm_free(hlu_data->gemm);
    SMPI_SHARED_FREE(hlu_data);
    SMPI_SHARED_FREE(platform_data);
end:
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

