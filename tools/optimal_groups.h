/*
 * find_optimal_group.h
 *
 *  Created on: 26 Mar 2015
 *      Author: xalid
 */

#ifndef TOOLS_OPTIMAL_GROUPS_H_
#define TOOLS_OPTIMAL_GROUPS_H_

#include "himpi/himpi.h"

int get_himpi_group(int msg_size, int root, MPI_Comm comm_world, int num_levels,
		int alg_in, int alg_out, himpi_operations op_id);

void save_himpi_optimal_groups(int min_msg_size, int max_msg_size,
		int msg_stride, int root, int num_levels, int alg_in, int alg_out,
		himpi_operations op_id, int use_one_proc, const char* file_name);

#endif /* TOOLS_OPTIMAL_GROUPS_H_ */
