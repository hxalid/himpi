/*
 * find_optimal_group.c
 *
 *  Created on: 26 Mar 2015
 *      Author: Khalid Hasanov
 */

#include "optimal_groups.h"
#include "MPIBlib/benchmarks/mpib.h"

#include <stdlib.h>

int get_hmpi_group(int count, MPI_Datatype datatype, int root,
		MPI_Comm comm_world, int num_levels, int alg_in, int alg_out,
		hmpi_operations op_id) {
	MPIB_result result;
	MPIB_precision precision;
	MPIB_getopt_precision_default(&precision);

	int my_rank;
	int num_procs;
	MPI_Comm_rank(comm_world, &my_rank);
	MPI_Comm_size(comm_world, &num_procs);

	int g, num_groups = 0;
	for (g = 1; g < num_procs; g++) {
		if (num_procs % g == 0)
			num_groups++;
	}

	double* g_times = (double*) calloc(num_groups, sizeof(double));
	if (g_times == NULL) {
		fprintf(stderr,
				"[get_hbcast_group]:Can't allocate memory for g_times\n");
		return -1; //TODO
	}

	MPI_Aint extent, lb;
	MPI_Type_get_extent(datatype, &lb, &extent);
	int message_size = extent * count;

	int i = 0;
	for (g = 1; g < num_procs; g++) {
		if (num_procs % g == 0) {
			MPIB_coll_container* container;
			if (op_id == op_bcast) {
				container = (MPIB_coll_container*) MPIB_HBcast_container_alloc(
						hierarchical_broadcast, g, num_levels, alg_in, alg_out);
			} else {
				container = (MPIB_coll_container*) MPIB_HReduce_container_alloc(
						hierarchical_reduce, g, num_levels, alg_in, alg_out);
			}

			int err = MPIB_measure_max(container, comm_world, 0, message_size,
					precision, &result);
			g_times[i++] = result.T;
		}
	}

	//TODO
	int group = gsl_stats_min_index(g_times, 1, num_groups);
	free(g_times);

	return group + 1;
}

/*
 * Calculate optimal number of groups for all number of processes
 * from HBCAST_MIN_PROCS up to comm_size and save it into a config file.
 */
void save_hmpi_optimal_groups(int count, MPI_Datatype datatype, int root,
		MPI_Comm comm_world, int num_levels, int alg_in, int alg_out, hmpi_operations op_id) {
	int g;
	int rank;
	int comm_size;
	int new_size;
	MPI_Comm_rank(comm_world, &rank);
	MPI_Comm_size(comm_world, &comm_size);

	FILE* fp;
	if (rank == 0) {
		//TODO: configurable filename
		fp = fopen(HMPI_CONF_FILE_NAME, "w"); //TODO:  should I overwrite?
		fprintf(fp, "#num_procs\tnum_groups\tnum_levels\talg_in\talg_out\n");
		if (fp == NULL) {
			printf("Try to open the configuration file %s\n", HMPI_CONF_FILE_NAME);
			perror("fopen");
			MPI_Abort(MPI_COMM_WORLD, 201);
		}
	}

	for (g = 0; g < comm_size - HBCAST_MIN_PROCS; g++) {
		MPI_Comm sub_comm;
		MPI_Comm_split(comm_world, (comm_size - rank > g) ? 0 : MPI_UNDEFINED,
				rank, &sub_comm);
		if (sub_comm != MPI_COMM_NULL) {
			MPI_Comm_size(sub_comm, &new_size);

			int group = get_hmpi_group(count, datatype, root, sub_comm,
					num_levels, alg_in, alg_out, op_id);
			if (group != -1) {
				if (rank == 0)
					fprintf(fp, "%d\t%d\t%d\t%d\t%d\n", new_size, group, 1, 0, 0);
			}

			MPI_Comm_free(&sub_comm);
		}
	}
	if (rank == 0)
		fclose(fp);

}

