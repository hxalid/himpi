/*
 * find_optimal_group.c
 *
 *  Created on: 26 Mar 2015
 *      Author: xalid
 */

#include "optimal_groups.h"

int get_hbcast_group(void *buffer, int count, MPI_Datatype datatype,
		int root, MPI_Comm comm_world, int num_groups, int rec, int alg){

	int group, res;

	res = MPI_HBcast(buffer, count, datatype, root, comm_world);

	return group;
}


