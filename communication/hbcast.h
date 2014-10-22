#ifndef HBCAST_H_
#define HBCAST_H_

#include <mpi.h>

#define MAX_INT_SIZE 50
#define HBCAST_MIN_PROCS 4
#define TEN_KB 10240
#define Bcast_TAG 904920477 
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


int validate_input(int num_groups, int num_procs);

t_bcast_response MPI_HBcast(void *buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm_world, int num_groups, int rec, int alg, int debug);

#endif
