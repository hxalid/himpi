/*
 * 2D LU factorization
 *
 * Authors: Khalid Hasanov
 */


#include "config.h"
#include "tools/hpnla_conf.h"
#include "tools/hpnla_debug.h"
#include "tools/hpnla_bind.h"
//#include "Matrix_init.h"
#include "lu_factorization.h"
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

    LU_data* lu_data = malloc(sizeof (LU_data));
    if (lu_data == NULL) {
        printf("malloc failed for lu_data\n");
        return -1;
    }
    
    Debug_data* debug_data = malloc(sizeof (Debug_data));
    if (debug_data == NULL) {
        printf("malloc failed for debug_data\n");
        return -1;
    }

    lu_data->block_size = 16;
    lu_data->m = 1024;
    lu_data->n = 1024;
    lu_data->k = 1024;
    lu_data->bcast_algorithm = 4; //original algorithm inside mpi
    lu_data->with_partial_pivoting = 1;

    size_t row, col;
    row = 0;
    col = 0;

    int opt, display = 0;
    char *conf_file = NULL;
    hpnla_process_conf conf;
    
    MPI_Init(&argc, &argv);

    lu_data->comm = MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &lu_data->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &lu_data->nproc);

    
    // for debug
#if DEBUG_MPI
    size_t loop = 1;
    while (loop == 1);
#endif

    optind = 1;
    
    //get the parameter from command line
    while ((opt = getopt(argc, argv, "hdsf:r:c:a:d:M:N:K:T:B:R:C:g:P:")) != -1) {
        switch (opt) {
            case 'h':
                if (lu_data->rank == 0)
                    info_print(0, row, col,
                        "Usage: LU_factorization [options]\n"
                        "	-M I	M size (default: %zu)\n"
                        "	-N I	N size (default: %zu)\n"
                        "	-K I	K size (default: %zu)\n"
                        "	-B I	Block size on the k dimension(default: %zu)\n"
                        "	-r I	processor row size (default: %zu)\n"
                        "	-c I	processor col size (default: %zu)\n"
                        "   -a I    broadcast algorithm (default: %d)\n"
                        "   -f {Filename} provide the file with the configuration\n"
                        "   -s  display the configuration file on the stderr\n"
                        "   -p I    should use partial pivoting (default: %d)\n"
                        "   -h	help\n",
                        lu_data->m,
                        lu_data->n,
                        lu_data->k,
                        lu_data->block_size,
                        row,
                        col,
                        lu_data->bcast_algorithm,
                        lu_data->with_partial_pivoting
                        );
                return 0;
            case 'M':
                lu_data->m = atoi(optarg);
                break;
            case 'N':
                lu_data->n = atoi(optarg);
                break;
            case 'K':
                lu_data->k = atoi(optarg);
                break;
            case 'B':
                lu_data->block_size = atoi(optarg);
                break;
            case 'r':
                lu_data->prow = atoi(optarg);
                break;
            case 'c':
                lu_data->pcol = atoi(optarg);
                break;
            case 'a':
                lu_data->bcast_algorithm = atoi(optarg);
                break;
            case 'p':
                lu_data->with_partial_pivoting = atoi(optarg);
                break;
            case 'f':
                conf_file = strdup(optarg);
                break;
            case 's':
                display = 1;
                break;
        }
    }
    
   // platform_data->nb_requested_procs = platform_data->size_row * platform_data->size_col;

    if (display == 1) {
        hpnla_print_conf(MPI_COMM_WORLD, lu_data->rank, stderr, "", "");
        MPI_Barrier(MPI_COMM_WORLD);
        exit(0);
    }

    if (conf_file == NULL) {
        info_print(0, row, col, "No configuration file\n");
    } else {
        conf = hpnla_get_conf(MPI_COMM_WORLD, conf_file);
    }

    if (conf.subopts[0] == '\0') {
        info_print(0, row, col, "No configuration for me inside the file\n");
    } else {
        info_print(0, row, col, "my_rank:%d, my_group:%s\n", conf.rank_intra, conf.subopts);
        if (hpnla_bind_process(conf.bind) == -1) goto end;
    }
    
    lu_data->gemm = hpnla_gemm_alloc(&conf);
    
    lu_factorize(lu_data, debug_data);

    // close the resources
    hpnla_gemm_free(lu_data->gemm);
    SMPI_SHARED_FREE(lu_data);
    SMPI_SHARED_FREE(debug_data);
end:
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    return 0;
}

