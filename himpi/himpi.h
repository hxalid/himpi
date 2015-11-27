/*!
 * himpi.h
 *
 *  Created on: 25 Mar 2015
 *      Author: Khalid Hasanov
 */

#ifndef HIMPI_HMPI_H_
#define HIMPI_HMPI_H_

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "config.h"

#define TEN_KB 10240
#define Bcast_TAG 904920477
#define ERR_GROUPS -99

#define HIMPI_ALG_IN 0
#define HIMPI_ALG_OUT 0
#define HIMPI_NUM_GROUPS 1
#define HIMPI_NUM_LEVELS 1
#define HIMPI_ROOT_PROC 0
#define HIMPI_MIN_PROCS 16
#define HIMPI_MIN_MSG 1024	   // 1kB
#define HIMPI_MAX_MSG 16777216  // 16MB
#define HIMPI_MSG_STRIDE 2
#define HIMPI_CONF_FILE "./default_file.conf"
#define HIMPI_GENERATE_CONFIG 0
#define _unused(x) ((void)x)

typedef enum himpi_operations {
	op_bcast, op_reduce, op_allreduce, op_scatter, op_gather, op_all
} himpi_operations;

/*! hmpi configuration parameters */
typedef struct hmpi_conf {
	int num_procs;
	int num_groups;
	int num_levels;
	int message_size;
	int alg_in;
	int alg_out;
	himpi_operations op_id;
} himpi_conf;

/*! himpi algorithm parameters */
typedef struct himpi_env {
	int min_msg;
	int max_msg;
	int msg_stride;
	int min_procs;
	int bcast_alg_in;
	int bcast_alg_out;
	int num_levels;
	int root;
	char* conf_file_name;
	int generate_config;
} himpi_env;

/*! It might be used in future to implement in-memory configuration */
typedef struct himpi_group_data {
	int num_procs;
	int num_groups;
	int msg_size;
} himpi_group_data;

extern const char *himpi_conf_file_name;
extern himpi_env henv;
extern himpi_group_data* group_data;
extern int himpi_debug;
extern int himpi_my_rank_world;
extern int himpi_num_ranks_world;
/* According to POSIX standard hostname should not exceed 255 bytes*/
extern char himpi_my_hostname[_POSIX_HOST_NAME_MAX];
extern MPI_Comm himpi_comm_world;

int validate_input(int num_groups, int num_procs);

int hierarchical_broadcast(void *buffer, int count, MPI_Datatype datatype,
		int root, MPI_Comm comm, int num_groups, int num_levels, int alg_in,
		int alg_out);

int hierarchical_reduce(void *sendbuf, void* recvbuf, int count,
		MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm_world,
		int num_groups, int num_levels, int alg_in, int alg_out);

int hierarchical_allreduce(void *sendbuf, void* recvbuf, int count,
		MPI_Datatype datatype, MPI_Op op, MPI_Comm comm_world, int num_groups,
		int num_levels, int alg_in, int alg_out);

int hierarchical_gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm, int num_groups, int num_levels, int alg_in, int alg_out);

int hierarchical_scatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm, int num_groups, int num_levels, int alg_in, int alg_out);

int HiMPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root,
		MPI_Comm comm_world);

int HiMPI_Reduce(void *sendbuf, void* recvbuf, int count,
		MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

int HiMPI_Allreduce(void *sendbuf, void* recvbuf, int count,
		MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int HiMPI_Scatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm);

int HiMPI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm);

int HiMPI_Init(int *argc, char ***argv);

#endif /* HIMPI_HMPI_H_ */
