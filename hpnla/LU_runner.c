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
#include "lu_factorization.h"
#include "hpnla_cblas.h"

#include <mpi.h>
#include <math.h>
#include <getopt.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char ** argv) {

    LU_data* lu_data= malloc(sizeof(LU_data));

    lu_data->Block_size_in = 32;
    lu_data->Block_size_out = lu_data->Block_size_in;
    lu_data->size_group_col = 1;
    lu_data->size_group_row = 1;
    lu_data->m = 1024;
    lu_data->n = 1024;
    lu_data->k = 1024;
    lu_data->bcast_algorithm = 4; //original algorithm inside mpi
    lu_data->distribution = 0; //default matrix distribution is non-cyclic

    int myrank;
    int NB_proc;
    size_t row, col;
    row = 0;
    col = 0;

    int opt, display = 0;
    char *conf_file = NULL;
    hpnla_process_conf conf;



    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &NB_proc);

    if (NB_proc != 1)
        for (lu_data->size_col = NB_proc / 2; NB_proc % lu_data->size_col; lu_data->size_col--);
    else
        lu_data->size_col = 1;

    lu_data->size_row = NB_proc / lu_data->size_col;
    if (lu_data->size_row > lu_data->size_col) {
        lu_data->size_col = lu_data->size_row;
        lu_data->size_row = NB_proc / lu_data->size_col;
    }


    // for the degub
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
                        lu_data->m,
                        lu_data->n,
                        lu_data->k,
                        lu_data->Block_size_in,
                        lu_data->Block_size_out,
                        lu_data->size_group_col,
                        lu_data->size_group_row,
                        lu_data->group,
                        lu_data->key,
                        row,
                        col,
                        lu_data->bcast_algorithm,
                        lu_data->distribution);
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
                lu_data->Block_size_in = atoi(optarg);
                break;
            case 'T':
                lu_data->Block_size_out = atoi(optarg);
                break;
            case 'R':
                lu_data->size_group_row = atoi(optarg);
                break;
            case 'C':
                lu_data->size_group_col = atoi(optarg);
                break;
            case 'g':
                lu_data->group = atoi(optarg);
                break;
            case 'k':
                lu_data->key = atoi(optarg);
                break;
            case 'r':
                lu_data->size_row = atoi(optarg);
                break;
            case 'c':
                lu_data->size_col = atoi(optarg);
                break;
            case 'a':
                lu_data->bcast_algorithm = atoi(optarg);
                break;
            case 'd':
                lu_data->distribution = atoi(optarg);
                break;
            case 'f':
                conf_file = strdup(optarg);
                break;
            case 's':
                display = 1;
                break;
        }
    }

    if (display == 1) {
        hpnla_print_conf(MPI_COMM_WORLD, myrank, stderr, "", "");
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

        info_print(0, row, col, "my_rank:%d, my_group:%s\n", conf.rank_intra, conf.subopts);

        if (hpnla_bind_process(conf.bind) == -1) goto end;


        lu_data->group = (size_t) atoi(conf.subopts);

    }

    lu_data->gemm = hpnla_gemm_alloc(&conf);

    lu_factorize(lu_data);

    // close properly the pragram
    hpnla_gemm_free(lu_data->gemm);

end:
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

