/*
 * hmpi.h
 *
 *  Created on: 25 Mar 2015
 *      Author: xalid
 */

#ifndef HMPI_HMPI_H_
#define HMPI_HMPI_H_

#include <mpi.h>
#include <stdio.h>

#include "config.h"

#define TEN_KB 10240
#define Bcast_TAG 904920477
#define ERR_GROUPS -99


#define HBCAST_MIN_PROCS 4
#define HREDUCE_MIN_PROCS 4


//extern int num_groups;
extern const char *HMPI_CONF_FILE_NAME;

int validate_input(int num_groups, int num_procs);

int HMPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm_world, int rec, int alg);


int HMPI_Reduce(void *snd_buffer, void* rcv_buffer, int count, MPI_Datatype datatype, MPI_Op op,
        int root, MPI_Comm comm_world, int num_groups, int rec, int alg, int debug, int myrank);

int HMPI_Init( int *argc, char ***argv );


#endif /* HMPI_HMPI_H_ */
