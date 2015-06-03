/*
 * find_optimal_group.h
 *
 *  Created on: 26 Mar 2015
 *      Author: xalid
 */

#ifndef TOOLS_OPTIMAL_GROUPS_H_
#define TOOLS_OPTIMAL_GROUPS_H_

#include "hmpi/hmpi.h"

int get_hmpi_group(int msg_size, int root, MPI_Comm comm_world, int num_levels,
		int alg_in, int alg_out, hmpi_operations op_id);

void save_hmpi_optimal_groups(int min_msg_size, int max_msg_size,
		int msg_stride, int root, MPI_Comm comm_world, int num_levels,
		int alg_in, int alg_out, hmpi_operations op_id, int use_one_proc);

#endif /* TOOLS_OPTIMAL_GROUPS_H_ */
