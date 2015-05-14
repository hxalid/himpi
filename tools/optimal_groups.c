/*
 * find_optimal_group.c
 *
 *  Created on: 26 Mar 2015
 *      Author: xalid
 */

#include "optimal_groups.h"
#include "MPIBlib/benchmarks/mpib.h"

#include <stdlib.h>

int get_hbcast_group(int count, MPI_Datatype datatype, int root,
		MPI_Comm comm_world, int rec, int alg) {

	MPIB_result result;
	MPIB_precision precision;
	MPIB_getopt_precision_default(&precision);
	MPIB_coll_container* container =
			(MPIB_coll_container*) MPIB_Bcast_container_alloc(HMPI_Bcast);

	int my_rank;
	int num_procs;
	MPI_Comm_rank(comm_world, &my_rank);
	MPI_Comm_size(comm_world, &num_procs);

	double* g_times = (double*) calloc(num_procs, sizeof(double));
	if (g_times == NULL) {
		fprintf(stderr,
				"[get_hbcast_group]:Can't allocate memory for g_times\n");
		return -1; //TODO
	}

	MPI_Aint extent, lb;
	MPI_Type_get_extent(datatype, &lb, &extent);
	int message_size = extent * count;

	int i = 0;
	for (i = 0; i < num_procs; i++) {
		int err = MPIB_measure_max(container, comm_world, 0, message_size,
				precision, &result);
		g_times[i] = result.T;
	}

	//TODO
	int group =
			my_rank == root ? gsl_stats_min_index(g_times, 1, num_procs) : -1;
	free(g_times);

	return group + 1; //TODO: num_groups = rank + 1
}

/*
 * Calculate optimal number of groups for all number of processes
 * from HBCAST_MIN_PROCS up to comm_size and save it into a config file.
 */
void save_hbcast_optimal_groups(int count, MPI_Datatype datatype, int root,
		MPI_Comm comm_world, int rec, int alg) {
	int i;
	int rank;
	int comm_size;
	int new_size;
	MPI_Comm_rank(comm_world, &rank);
	MPI_Comm_size(comm_world, &comm_size);

	//TODO
	for (i = 1; i < comm_size-HBCAST_MIN_PROCS; i++) {
		MPI_Comm sub_comm;
		MPI_Comm_split(comm_world, (comm_size - rank > i) ? 0 : MPI_UNDEFINED,
				rank, &sub_comm);
		if (sub_comm != MPI_COMM_NULL) {
			MPI_Comm_size(sub_comm, &new_size);

			int group = get_hbcast_group(count, datatype, root, sub_comm, rec,
					alg);
			if (group != -1) {
				//TODO: write (comm_size, new_size, group) into a config file
				fprintf(stdout, "comm_size=%d, rank=%d, new_size=%d, i=%d, group=%d\n",
						comm_size, rank, new_size, i, group);
			}

			MPI_Comm_free(&sub_comm);
		}
	}

}

