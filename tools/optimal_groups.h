/*
 * find_optimal_group.h
 *
 *  Created on: 26 Mar 2015
 *      Author: xalid
 */

#ifndef TOOLS_OPTIMAL_GROUPS_H_
#define TOOLS_OPTIMAL_GROUPS_H_

#include "hmpi/hmpi.h"

typedef enum hmpi_operations {
   op_bcast,
   op_reduce
} hmpi_operations;

int get_hmpi_group(int count, MPI_Datatype datatype, int root,
		MPI_Comm comm_world, int num_levels, int alg_in, int alg_out, hmpi_operations op_id);

void save_hbcast_optimal_groups(int count, MPI_Datatype datatype, int root,
		MPI_Comm comm_world, int num_levels, int alg_in, int alg_out, hmpi_operations op_id);

#endif /* TOOLS_OPTIMAL_GROUPS_H_ */
