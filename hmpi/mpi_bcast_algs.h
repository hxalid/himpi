#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_INT_SIZE 50
#define TEN_KB 10240
#define Bcast_TAG 904920477 
#define bcast_linear_segment_size 1024

typedef enum hpnla_bcast_algo {
    lin,
    binary, // 1
    flat,
    binomial,
    original, //4
    scatter_lr_allgather,
    scatter_rd_allgather,
    binomial_mpich, //7,
    pipelined_linear // 8
} hpnla_bcast_algo;




void bcast_scatter_lr_allgather(void * buff, int count, MPI_Datatype data_type,
        int root, MPI_Comm comm);

int bcast_scatter_rdb_allgather(void * buff, int count, MPI_Datatype data_type,
        int root, MPI_Comm comm);

int bcast_binomial_tree(void * buff, int count, MPI_Datatype data_type,
        int root, MPI_Comm comm);

int bcast_pipelined_linear(void *buf, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm);
