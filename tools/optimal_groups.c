/*
 * find_optimal_group.c
 *
 *  Created on: 26 Mar 2015
 *      Author: Khalid Hasanov
 */

#include "optimal_groups.h"
#include "MPIBlib/benchmarks/mpib.h"

#include <stdlib.h>

hmpi_group_data* group_data;

int get_hmpi_group(int msg_size, int root,
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

			int err = MPIB_measure_max(container, comm_world, 0, msg_size,
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
void save_hmpi_optimal_groups(int min_msg_size, int max_msg_size, int msg_stride, int root,
		MPI_Comm comm_world, int num_levels, int alg_in, int alg_out, hmpi_operations op_id, int use_one_proc) {
	int p;
	int rank;
	int comm_size;
	int new_size;
	MPI_Comm_rank(comm_world, &rank);
	MPI_Comm_size(comm_world, &comm_size);

	FILE* fp;
	if (rank == 0) {
		//TODO: configurable filename
		fp = fopen(HMPI_CONF_FILE_NAME, "w"); //TODO:  should I overwrite?
		fprintf(fp, "#num_procs\tnum_groups\tnum_levels\tmsg_size\talg_in\talg_out\n");
		if (fp == NULL) {
			fprintf(stdout, "Try to open the configuration file %s\n", HMPI_CONF_FILE_NAME);
			perror("fopen");
			MPI_Abort(MPI_COMM_WORLD, 201);
		}
	}

	/*
	 * Find optimal number of groups for p\in[HMPI_MIN_PROCS+1, comm_size]
	 */
	int p_end = comm_size - HMPI_MIN_PROCS;
	if (use_one_proc) {
		p_end = 1;
	}

	for (p = 0; p < p_end; p++) {
		MPI_Comm sub_comm;
		MPI_Comm_split(comm_world, (comm_size - rank > p) ? 0 : MPI_UNDEFINED,
				rank, &sub_comm);
		if (sub_comm != MPI_COMM_NULL) {
			MPI_Comm_size(sub_comm, &new_size);

			int msg = 0;
			for(msg = min_msg_size; msg<=max_msg_size; msg*=msg_stride) {

				int group = get_hmpi_group(msg, root, sub_comm,
						num_levels, alg_in, alg_out, op_id);
				if (group != -1) {
					if (rank == 0)
						fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\n", new_size, group, 1, msg, 0, 0);
				}
			}

			MPI_Comm_free(&sub_comm);
		}
	}

	if (rank == 0)
		fclose(fp);

}


/*
 * TODO: not ready yet
 */
void save_hmpi_groups_in_memory(int min_msg_size, int max_msg_size, int msg_stride, int root,
		MPI_Comm comm_world, int num_levels, int alg_in, int alg_out, hmpi_operations op_id, int use_one_proc) {
	int p;
	int rank;
	int comm_size;
	int new_size;
	MPI_Comm_rank(comm_world, &rank);
	MPI_Comm_size(comm_world, &comm_size);

	int p_end = comm_size - HMPI_MIN_PROCS;
	if (use_one_proc) {
		p_end = 1;
	}

	int mc = 0;
	for(mc = min_msg_size; mc<=max_msg_size; mc*=msg_stride) {
		mc++;
	}

	int data_size = p_end * mc;
	group_data = (hmpi_group_data*)malloc(data_size*sizeof(hmpi_group_data));

	int i = 0;

	for (p = 0; p < p_end; p++) {
			MPI_Comm sub_comm;
			MPI_Comm_split(comm_world, (comm_size - rank > p) ? 0 : MPI_UNDEFINED,
					rank, &sub_comm);
			if (sub_comm != MPI_COMM_NULL) {
				MPI_Comm_size(sub_comm, &new_size);

				int msg = 0;
				for(msg = min_msg_size; msg<=max_msg_size; msg*=msg_stride) {

					int group = get_hmpi_group(msg, root, sub_comm,
							num_levels, alg_in, alg_out, op_id);
					if (group != -1) {
						group_data[i].num_procs = new_size;
						group_data[i].msg_size = msg;
						group_data[i].num_groups = group;
					}
				}

				MPI_Comm_free(&sub_comm);
			}
	}
}


int group_data_cmp (const struct hmpi_group_data *gd1, const struct hmpi_group_data *gd2){
  return (gd1->num_procs==gd2->num_procs && gd1->num_groups==gd2->num_groups && gd1->msg_size==gd2->msg_size);
}

hmpi_group_data* find_group_data(int num_procs, int num_groups, int msg_size) {
	struct hmpi_group_data needed_data, *result;
	needed_data.num_procs = num_procs;
	needed_data.num_groups = num_groups;
	needed_data.msg_size = msg_size;

	int count = sizeof(group_data)/sizeof(struct hmpi_group_data);

	result = bsearch (&needed_data, group_data, count, sizeof (struct hmpi_group_data),
			group_data_cmp);

	return result;
}











