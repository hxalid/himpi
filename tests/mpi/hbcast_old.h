#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define HBCAST_MIN_PROCS 4
#define ERR_GROUPS -99


typedef struct {
    int rec_out_group_comm;
    int rec_in_group_comm;
    int rec_comm_world;
    int err_flag;
    double in_time;
    double out_time;
    double group_time;
    double no_group_time;
} t_bcast_response;


t_bcast_response MPI_HBcast_old(void *buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm_world, int num_groups, int rec, int alg, int debug);


t_bcast_response MPI_HBcast_Opt(void *buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm_world, MPI_Comm in_group_comm,
        MPI_Comm out_group_comm, int num_groups, int rec, int alg, int debug);
