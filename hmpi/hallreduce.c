/*
 * hallreduce.c
 *
 *  Created on: 26 Mar 2015
 *      Author: Khalid Hasanov
 */

#include "hmpi.h"

int HMPI_Allreduce(void *snd_buffer, void* rcv_buffer, int count,
		MPI_Datatype datatype, MPI_Op op, MPI_Comm comm_world, int num_groups,
		int rec, int alg, int debug, int myrank) {
	int pg;
	int rank;
	int size;
	int stride;
	int root_inside;
	int root_outside;
	int my_group;
	double *reduce_in;

	MPI_Comm in_group_comm, out_group_comm;

	MPI_Comm_rank(comm_world, &rank);
	MPI_Comm_size(comm_world, &size);
	if (size == 1)
		return MPI_SUCCESS;

	/*TODO make num_groups configurable*/
	if (size > HREDUCE_MIN_PROCS && (num_groups > 1 && num_groups < size)) {
		pg = size / num_groups;
		my_group = rank / pg;

		MPI_Comm_split(comm_world, my_group, rank, &in_group_comm);
		// MPI_Comm_split(comm_world, (rank - my_group * pg == 0) ? 0 : MPI_UNDEFINED, rank, &out_group_comm);
		MPI_Comm_split(comm_world, rank - my_group * pg, rank, &out_group_comm);

		if (rec == 1)
			reduce_in = (double*) malloc(count * sizeof(double));

		//Start broadcast inside groups
		switch (rec) {
		case 1: // 1 level of hierarchy
			MPI_Allreduce(snd_buffer, reduce_in, count, datatype, op, in_group_comm);
			break;
		case -1:
			// Just to see if broadcast inside groups is better than that of with MPI_COMM_WORLD
			MPI_Allreduce(snd_buffer, rcv_buffer, count, datatype, op, in_group_comm);
			break;
		default: // e.g. -2
			// Don't broadcast inside groups. Just to see if broadcast between groups is better than that of with MPI_COMM_WORLD
			break;
		}

		//Start allreduce between groups
		if (out_group_comm != MPI_COMM_NULL && rec != -1) { // if rec == -1 then allreduce only inside
			int out_size = 0;
			int out_rank = -1;

			MPI_Allreduce((rec == 1) ? reduce_in : snd_buffer, rcv_buffer, count, datatype, op, out_group_comm);
		}

		// Now broadcast sub-sums inside groups. We need it if we use communicators between the groups leaders
		/*  if (rec == 1) {
		 MPI_Bcast(rcv_buffer, count, datatype, 0, in_group_comm);
		 }
		 */

		//TODO: check
		if (rec == 1)
			free(reduce_in);

		if (out_group_comm != MPI_COMM_NULL)
			MPI_Comm_free(&out_group_comm);
		if (in_group_comm != MPI_COMM_NULL)
			MPI_Comm_free(&in_group_comm);
	} else if (size <= HREDUCE_MIN_PROCS || (num_groups == 1 || num_groups == size)) {
		//  fprintf(stdout, "Using non-hierarchical reduce: [p=%d, g=%d]\n", size, num_groups);
		MPI_Allreduce(snd_buffer, rcv_buffer, count, datatype, op, comm_world);
	} else {
		/*TODO*/
		// fprintf(stdout, "Wrong number of groups: [p=%d, g=%d]\n", size, num_groups);
	}

	return MPI_SUCCESS;
}

