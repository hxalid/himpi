/*!
 * hmpi.h
 *
 *  Created on: 25 Mar 2015
 *      Author: Khalid Hasanov
 */

#ifndef HMPI_HMPI_H_
#define HMPI_HMPI_H_

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "config.h"

#define TEN_KB 10240
#define Bcast_TAG 904920477
#define ERR_GROUPS -99

#define HMPI_ALG_IN 0
#define HMPI_ALG_OUT 0
#define HMPI_NUM_GROUPS 1
#define HMPI_NUM_LEVELS 1
#define HMPI_ROOT_PROC 0
#define HMPI_MIN_PROCS 4
#define HMPI_MIN_MSG 1024	   // 1Kb
#define HMPI_MAX_MSG 16777216  // 16Mb
#define HMPI_MSG_STRIDE 2
#define HMPI_CONF_FILE "./default_file.conf"
#define HMPI_GENERATE_CONFIG 0
#define _unused(x) ((void)x)

typedef enum hmpi_operations {
	op_bcast, op_reduce, op_allreduce, op_scatter, op_gather
} hmpi_operations;

/*! hmpi configuration parameters */
typedef struct hmpi_conf {
	int num_procs;
	int num_groups;
	int num_levels;
	int message_size;
	int alg_in;
	int alg_out;
} hmpi_conf;

/*! hmpi algorithm parameters */
typedef struct hmpi_env {
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
} hmpi_env;

/*! It might be used in future to implement in-memory configuration */
typedef struct hmpi_group_data {
	int num_procs;
	int num_groups;
	int msg_size;
} hmpi_group_data;

extern const char *HMPI_CONF_FILE_NAME;
extern hmpi_env henv;
extern hmpi_group_data* group_data;

int validate_input(int num_groups, int num_procs);

int hierarchical_broadcast(void *buffer, int count, MPI_Datatype datatype,
		int root, MPI_Comm comm, int num_groups, int num_levels, int alg_in,
		int alg_out);

int hierarchical_reduce(void *snd_buffer, void* rcv_buffer, int count,
		MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm_world,
		int num_groups, int num_levels, int alg_in, int alg_out);

int hierarchical_allreduce(void *snd_buffer, void* rcv_buffer, int count,
		MPI_Datatype datatype, MPI_Op op, MPI_Comm comm_world, int num_groups,
		int num_levels, int alg_in, int alg_out);

int hierarchical_gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm, int num_groups, int num_levels, int alg_in, int alg_out);

int hierarchical_scatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm, int num_groups, int num_levels, int alg_in, int alg_out);

int HMPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root,
		MPI_Comm comm_world);

int HMPI_Reduce(void *snd_buffer, void* rcv_buffer, int count,
		MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

int HMPI_Allreduce(void *snd_buffer, void* rcv_buffer, int count,
		MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int HMPI_Scatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm);

int HMPI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm);

int HMPI_Init(int *argc, char ***argv);

#endif /* HMPI_HMPI_H_ */
