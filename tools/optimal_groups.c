/*
 * find_optimal_group.c
 *
 *  Created on: 26 Mar 2015
 *      Author: xalid
 */

#include "optimal_groups.h"
#include "MPIBlib/benchmarks/mpib.h"

#include <stdlib.h>

int get_hbcast_group(int count, MPI_Datatype datatype,
		int root, MPI_Comm comm_world, int num_groups, int rec, int alg){

	MPIB_result result;
	MPIB_precision precision;
	MPIB_getopt_precision_default(&precision);
	MPIB_coll_container* container = (MPIB_coll_container*)MPIB_Bcast_container_alloc(HMPI_Bcast);

	int my_rank;
	MPI_Comm_rank(comm_world, &my_rank);

	double* g_times = (double*)calloc(num_groups, sizeof(double));
	if (g_times == NULL) {
		fprintf(stderr, "[get_hbcast_group]:Can't allocate memory for g_times\n");
		return -1; //TODO
	}

	MPI_Aint extent, lb;
	MPI_Type_get_extent(datatype, &lb, &extent);
	int message_size =extent*count;

	int i=0;
	for(i=0; i<num_groups; i++) {
		int err = MPIB_measure_max(container, comm_world, 0, message_size, precision, &result);
		g_times[i] = result.T;
	}

	//TODO
	int group = my_rank==root?gsl_stats_min_index(g_times, 1, num_groups):-1;
	free(g_times);

	return group + 1; //TODO: num_groups = rank + 1
}


