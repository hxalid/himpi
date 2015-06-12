/*!
 * hscatter.c
 *
 *  Created on: 25 Mar 2015
 *      Author: Khalid Hasanov
 */

#include "hmpi.h"
#include "tools/utils.h"


int hierarchical_scatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm, int num_groups, int num_levels, int alg_in, int alg_out) {
	_unused(num_levels);
	_unused(alg_in );
	_unused(alg_out);

	int res = MPI_SUCCESS;

	if (num_groups == 1) {
		res = MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt,
				recvtype, root, comm);
	} else {
		int num_procs;
		int rank;
		MPI_Comm in_group_comm, out_group_comm;

		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &num_procs);

		if (num_procs % num_groups != 0 || num_procs <= num_groups) {
			fprintf(stderr, "HScatter: Wrong number of groups\n");
			return -1;
		}

		int pg = num_procs / num_groups;
		int my_group = rank / pg;

		MPI_Aint extent, lb;
		MPI_Type_get_extent(recvtype, &lb, &extent);

		void* g_rcvbuf = malloc(extent * pg * recvcnt);

		/*! Scatter between groups */
		MPI_Comm_split(comm, (rank - my_group * pg == 0) ? 0 : MPI_UNDEFINED,
				rank, &out_group_comm);
		if (out_group_comm != MPI_COMM_NULL) {
			res = MPI_Scatter(sendbuf, pg * recvcnt, sendtype, g_rcvbuf,
					pg * recvcnt, recvtype, root, out_group_comm);
			MPI_Comm_free(&out_group_comm);
		}

		//TODO
		if (res != MPI_SUCCESS) {
			MPI_Abort(comm, 1);
		}

		/*! Scatter inside groups */
		MPI_Comm_split(comm, my_group, rank, &in_group_comm);
		if (in_group_comm != MPI_COMM_NULL) {
			res = MPI_Scatter(g_rcvbuf, recvcnt, sendtype, recvbuf, recvcnt,
					recvtype, root, in_group_comm);
			MPI_Comm_free(&in_group_comm);
		}

		free(g_rcvbuf);
	}
	return res;
}

int HMPI_Scatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm) {
	MPI_Aint extent, lb;
	MPI_Type_get_extent(sendtype, &lb, &extent);
	int msg_size = extent * sendcnt;

	hmpi_conf my_conf = hmpi_get_my_conf(comm, msg_size, root,
			hmpi_conf_file_name, op_scatter);

	return hierarchical_scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt,
			recvtype, root, comm, my_conf.num_groups, my_conf.num_levels, my_conf.alg_in, my_conf.alg_out);
}










