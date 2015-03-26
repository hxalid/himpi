/*
 * hscatter.c
 *
 *  Created on: 25 Mar 2015
 *      Author: Khalid Hasanov
 */

#include "hmpi.h"

/*
 make the send and receive buffer types general
 */
int MPI_HScatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm, int num_groups) {

	if (num_groups == 1) {
		return MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt,
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

//		int data_type_size;
//		MPI_Type_size(recvtype, &data_type_size);
		MPI_Aint extent, lb;
		MPI_Type_get_extent(recvtype, &lb, &extent);

//		void* g_rcvbuf = malloc(data_type_size * pg * recvcnt);
		void* g_rcvbuf = malloc(extent * pg * recvcnt);

		MPI_Comm_split(comm, (rank - my_group * pg == 0) ? 0 : MPI_UNDEFINED,
				rank, &out_group_comm);
		MPI_Comm_split(comm, my_group, rank, &in_group_comm);

		/* Scatter between groups */
		if (out_group_comm != MPI_COMM_NULL) {
			MPI_Scatter(sendbuf, pg * recvcnt, sendtype, g_rcvbuf, pg * recvcnt,
					recvtype, root, out_group_comm);
			MPI_Comm_free(&out_group_comm);
		}

		/* Scatter inside groups */
		if (in_group_comm != MPI_COMM_NULL) {
			MPI_Scatter(g_rcvbuf, recvcnt, sendtype, recvbuf, recvcnt, recvtype,
					root, in_group_comm);
			MPI_Comm_free(&in_group_comm);
		}

		free(g_rcvbuf);

		return MPI_SUCCESS;
	}
}

