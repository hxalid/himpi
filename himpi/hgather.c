/*
 * hgather.c
 *
 *  Created on: 26 Mar 2015
 *      Author: Khalid Hasanov
 */

#include "himpi.h"
#include "tools/utils.h"


int hierarchical_gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm, int num_groups, int num_levels, int alg_in, int alg_out) {
	_unused(num_levels);
	_unused(alg_in);
	_unused(alg_out);

	int res = MPI_SUCCESS;

	if (num_groups == 1) {
		res = MPI_Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype,
				root, comm);
	} else {
		int rank;
		int num_procs;
		MPI_Comm in_group_comm, out_group_comm;

		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &num_procs);

		int pg = num_procs / num_groups;
		int my_group = rank / pg;

		if (num_procs % num_groups != 0 || num_procs <= num_groups) {
			fprintf(stderr, "MPI_HGather: Wrong number of groups\n");
			return -1;
		}

		MPI_Aint extent, lb;
		MPI_Type_get_extent(recvtype, &lb, &extent);
		void* g_rcvbuf = malloc(extent * pg * recvcnt);

		/*! Gather inside groups */
		MPI_Comm_split(comm, my_group, rank, &in_group_comm);
		if (in_group_comm != MPI_COMM_NULL) {
			res = MPI_Gather(sendbuf, recvcnt, sendtype, g_rcvbuf, recvcnt,
					sendtype, root, in_group_comm);
			MPI_Comm_free(&in_group_comm);
		}

		//TODO
		if (res != MPI_SUCCESS) {
			MPI_Abort(comm, 1);
		}

		/*! Gather between groups */
		MPI_Comm_split(comm, (rank - my_group * pg == 0) ? 0 : MPI_UNDEFINED,
				rank, &out_group_comm);
		if (out_group_comm != MPI_COMM_NULL) {
			res = MPI_Gather(g_rcvbuf, pg * recvcnt, sendtype, recvbuf,
					pg * recvcnt, sendtype, root, out_group_comm);
			MPI_Comm_free(&out_group_comm);
		}

		free(g_rcvbuf);
	}

	return res;
}

int HiMPI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm) {
	MPI_Aint extent, lb;
	MPI_Type_get_extent(sendtype, &lb, &extent);
	int msg_size = extent * sendcnt;

	himpi_conf my_conf = himpi_get_my_conf(comm, msg_size, root,
			himpi_conf_file_name, op_gather);

	return hierarchical_gather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt,
			recvtype, root, comm, my_conf.num_groups, my_conf.num_levels,
			my_conf.alg_in, my_conf.alg_out);
}

