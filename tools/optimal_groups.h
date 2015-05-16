/*
 * find_optimal_group.h
 *
 *  Created on: 26 Mar 2015
 *      Author: xalid
 */

#ifndef TOOLS_OPTIMAL_GROUPS_H_
#define TOOLS_OPTIMAL_GROUPS_H_

#include "hmpi/hmpi.h"

int get_hbcast_group(int count, MPI_Datatype datatype, int root,
		MPI_Comm comm_world, int num_levels, int alg_in, int alg_out);

void save_hbcast_optimal_groups(int count, MPI_Datatype datatype, int root,
		MPI_Comm comm_world, int num_levels, int alg_in, int alg_out);

#endif /* TOOLS_OPTIMAL_GROUPS_H_ */
