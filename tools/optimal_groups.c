/*!
 * find_optimal_group.c
 *
 *  Created on: 26 Mar 2015
 *      Author: Khalid Hasanov
 */

#include "optimal_groups.h"
#include "utils.h"
#include "MPIBlib/benchmarks/mpib.h"

#include <unistd.h>
#include <stdlib.h>
#include <math.h>

himpi_group_data* group_data;

//TODO:  define new config option to set if we have parallel file system
//#define HIMPI_HAS_PFS 1

/*!
 *  This function is used to find the requested
 *  factor of the given process/number.
 */
int get_specific_factor(int num_procs, int factor_idx) {
	int g = 0;
	int idx_counter = 0;
	int factor = 1;

	if (factor_idx > 0) {
		for (g = 1; g < num_procs; g++) {
			if (num_procs % g == 0) {
				idx_counter++;
				if (idx_counter == factor_idx) {
					factor = g;
					break;
				}
			}
		}
	}
	return factor;
}

/*!
 * \param msg_size a message size
 * \param root rank of the root process
 * \param comm_world parent communicator from which all sub-communicators will be created
 * \param num_levels number of hierarchies
 * \param alg_in collective algorithm id to use inside groups
 * \param alg_out collective algorithms id to use among groups
 * \param op_id hmpi collective operation id to use. IDs are defined in hmpi_operation in hmpi.h
 * \return optimal number of groups.
 */
int get_himpi_group(int msg_size, int root, MPI_Comm comm_world, int num_levels,
		int alg_in, int alg_out, himpi_operations op_id) {
	MPIB_result result;
	MPIB_precision precision;
	MPIB_getopt_precision_default(&precision);

	int num_procs;
	MPI_Comm_size(comm_world, &num_procs);

	int g, num_groups = 0;
	for (g = 1; g < num_procs; g++) {
		if (num_procs % g == 0)
			num_groups++;
	}

	double* g_times = (double*) calloc(num_groups, sizeof(double));
	if (g_times == NULL) {
		himpi_err("[get_hbcast_group]:Can't allocate memory for g_times\n");
		return -1; //TODO
	}

	int i = 0;
	for (g = 1; g < num_procs; g++) {
		if (num_procs % g == 0) {
			MPIB_coll_container* container;
			switch (op_id) {
			case op_bcast:
				container = (MPIB_coll_container*) MPIB_HBcast_container_alloc(
						hierarchical_broadcast, g, num_levels, alg_in, alg_out);
				break;
			case op_reduce:
				container = (MPIB_coll_container*) MPIB_HReduce_container_alloc(
						hierarchical_reduce, g, num_levels, alg_in, alg_out);
				break;
			case op_allreduce:
				container =
						(MPIB_coll_container*) MPIB_HAllreduce_container_alloc(
								hierarchical_allreduce, g, num_levels, alg_in,
								alg_out);
				break;
			case op_scatter:
				container =
						(MPIB_coll_container*) MPIB_HScatter_container_alloc(
								hierarchical_scatter, g, num_levels, alg_in,
								alg_out);
				break;
			case op_gather:
				container = (MPIB_coll_container*) MPIB_HGather_container_alloc(
						hierarchical_gather, g, num_levels, alg_in, alg_out);
				break;
			default:
				fprintf(stdout, "Unknown operation id: %d\n", op_id);
				MPI_Abort(MPI_COMM_WORLD, 201);  //TODO
				break;
			}

			int err = MPIB_measure_max(container, comm_world, root, msg_size,
					precision, &result);
			g_times[i++] = result.T;
		}
	}

	int min_idx = gsl_stats_min_index(g_times, 1, num_groups);
	/*!
	 * Send min_idx+1 as array index start from 0,
	 * but we prefer to start from 1 to index factors.
	 */
	int group = get_specific_factor(num_procs, min_idx + 1);
	free(g_times);

	return group;
}

/*!
 * Calculate optimal number of groups for all number of processes
 * from HBCAST_MIN_PROCS up to comm_size and save it into a config file.
 */
void save_himpi_optimal_groups(int min_msg_size, int max_msg_size,
		int msg_stride, int root, int num_levels, int alg_in, int alg_out,
		himpi_operations op_id, int use_one_proc, const char* file_name) {
	int p;
	int new_size;
	int file_exists = 0;

	/*!
	 * Find optimal number of groups for p\in[HMPI_MIN_PROCS+1, himpi_num_ranks_world]
	 */
	int p_start = HIMPI_MIN_PROCS;
	if (use_one_proc) {
		p_start = himpi_num_ranks_world;
	}

	char* config_file_name = create_file_name(file_name, op_id);

	if (is_same_config(min_msg_size, max_msg_size, msg_stride,
			himpi_num_ranks_world, op_id, config_file_name))
		return;

	FILE* fp;
	if( access( config_file_name, F_OK ) != -1 ) {
	   fp = fopen(config_file_name, "a"); //File exists
	   file_exists = 1;
	} else
	   fp = fopen(config_file_name, "w");  // File does not exist

	//free(config_file_name);

	if (fp == NULL) {
		himpi_abort(-1,
				"save_himpi_optimal_groups: Can't open configuration file: %s @ %s:%d",
				himpi_conf_file_name, __FILE__, __LINE__);
	}

	/*
	 *  If we have a parallel file system then only one rank writes to file,
	 *  otherwise all create there own files.
	 */
	if ( !file_exists )
	#ifdef HIMPI_HAS_PFS
	    if ( !himpi_my_rank_world )
	       fprintf(fp, "#num_procs\tnum_groups\tnum_levels\tmsg_size\talg_in\talg_out\top_id\n");
    #else
	    fprintf(fp, "#num_procs\tnum_groups\tnum_levels\tmsg_size\talg_in\talg_out\top_id\n");
    #endif

	for (p = himpi_num_ranks_world; p >= p_start; p /= 2) {
		//TODO:  do we need all number of processes? May be we can use only pow of 2 number of processes.
		MPI_Comm sub_comm;
		MPI_Comm_split(himpi_comm_world,
				(himpi_my_rank_world < p) ? 0 : MPI_UNDEFINED,
				himpi_my_rank_world, &sub_comm);
		if (sub_comm != MPI_COMM_NULL) {
			MPI_Comm_size(sub_comm, &new_size);

			int msg = 0;
			for (msg = min_msg_size; msg <= max_msg_size; msg *= msg_stride) {
				int group = get_himpi_group(msg, root, sub_comm, num_levels,
						alg_in, alg_out, op_id);

				if (group != -1) {
					#if HIMPI_HAS_PFS
						if ( !himpi_my_rank_world ) {
							fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\n", new_size, group, 1, msg, 0, 0, op_id);
							printf(" %d\t%d\t%d\t%d\t%d\t%d\t%d\n", new_size, group, 1, msg, 0, 0, op_id);
						}
					#else
					    fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\n", new_size, group, 1, msg, 0, 0, op_id);
					#endif

					fflush(fp);
				}
			}

			MPI_Comm_free(&sub_comm);
		}
	}

	/*
	 *  TODO:  think twice if we need barrier before closing file
	 */
	fclose(fp);
}

/*!
 * TODO: not implemented yet
 */
void save_hmpi_groups_in_memory(int min_msg_size, int max_msg_size,
		int msg_stride, int root, int num_levels, int alg_in, int alg_out,
		himpi_operations op_id, int use_one_proc) {
	int p;
	int comm_size;
	int new_size;

	int p_end = comm_size - HIMPI_MIN_PROCS;
	if (use_one_proc) {
		p_end = 1;
	}

	int mc = 0;
	for (mc = min_msg_size; mc <= max_msg_size; mc *= msg_stride) {
		mc++;
	}

	int data_size = p_end * mc;
	group_data = (himpi_group_data*) malloc(
			data_size * sizeof(himpi_group_data));

	int i = 0;

	for (p = 0; p < p_end; p++) {
		MPI_Comm sub_comm;
		MPI_Comm_split(himpi_comm_world,
				(comm_size - himpi_my_rank_world > p) ? 0 : MPI_UNDEFINED,
				himpi_my_rank_world, &sub_comm);
		if (sub_comm != MPI_COMM_NULL) {
			MPI_Comm_size(sub_comm, &new_size);

			int msg = 0;
			for (msg = min_msg_size; msg <= max_msg_size; msg *= msg_stride) {

				int group = get_himpi_group(msg, root, sub_comm, num_levels,
						alg_in, alg_out, op_id);
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

/*
 int group_data_cmp(const struct hmpi_group_data *gd1,
 const struct hmpi_group_data *gd2) {
 return (gd1->num_procs == gd2->num_procs
 && gd1->num_groups == gd2->num_groups
 && gd1->msg_size == gd2->msg_size);
 }

 hmpi_group_data* find_group_data(int num_procs, int num_groups, int msg_size) {
 struct hmpi_group_data needed_data, *result;
 needed_data.num_procs = num_procs;
 needed_data.num_groups = num_groups;
 needed_data.msg_size = msg_size;

 int count = sizeof(group_data) / sizeof(struct hmpi_group_data);

 result = bsearch(&needed_data, group_data, count,
 sizeof(struct hmpi_group_data), group_data_cmp);

 return result;
 }
 */

