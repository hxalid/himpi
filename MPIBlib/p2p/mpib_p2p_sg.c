#include "mpib_p2p.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

extern char *all_names;
extern int original_comm_size;
extern int PROC_NUM_SPAWNED;
extern MPI_Comm intracomm;
extern MPI_Comm control_intracomm;

extern MPI_Request* global_reqs;
extern int global_reqs_counter;

/*!
 * Send out a message to all spawned MPI processes to let them know
 * how much data to forward:
 * if flag = 1: then this is a modified point-to-point communication and
 * 	the rank_list contains all processes on the two communicating nodes
 * if flag = 2: then this is a terminating call to all slave processes on all nodes
 * The call is blocking since it is usually followed by non-blocking communication which depends on this synchronisation
 */
int signal_spawned_procs(/*int* rank_list,*/ int size, int flag, int source, int target, int sendcount, int rest) {
	MPI_Request reqs[size-original_comm_size];
	int info[4];
	info[0] = flag;
	info[1] = target;
	info[2] = sendcount;
	info[3] = rest;
	int i;
	for (i = original_comm_size; i < size; i++) {
			MPI_Isend(&info, 4, MPI_INT, i, 0, control_intracomm,
					&reqs[i-original_comm_size]);
	}
	MPI_Waitall(size-original_comm_size, reqs, MPI_STATUSES_IGNORE);
	return 0;
}

/*!
 * This method is the scatter-gather implementation of point-to-point communication. It should be called by both sender and receiver -
 * other than actual MPI_Send/MPI_Recv calls. This means there are strictly source/target/rest processes related sections. Also, this method
 * can be both non-blocking or blocking. The non-blocking version keeps an internal global request list and must be used in combination with MPIB_Waitall. The blocking version does not and waits for completion.
 *
 * Source process:
 * - tells the spawned processes at the sender/receiver nodes to participate in blocking fashion
 * - it sends around data in a non-blocking way; if one data chunk is of irregular size, it goes to the target process
 *
 * Target process:
 * - blockingly receives a data chunk (potentially of irregular size) from the source process
 * - non-blockingly gets chunks from everyone
 *
 * Spawned processes: blockingly receive and send a chunk of data
 */
int MPIB_scatter_gather_based_p2p(void* sendbuf, int sendcount, int rest,
		MPI_Datatype sendtype, void* recvbuf, int recvcount,
		MPI_Datatype recvtype,  MPI_Comm comm, int source, int target, void* interbuf, int MPIB_p2p_blocking) {

	int rank, size;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);
	if (rank == source) {
		// The source node signals the spawned processes
		// on the sender and receiver node to participate
		signal_spawned_procs(size, 1, source, target, sendcount, rest);
		// Allocate the requests both for scatter and gather operation
		//MPI_Request reqs[size-original_comm_size+1];
		global_reqs = (MPI_Request *) realloc(global_reqs,
				(global_reqs_counter+size-original_comm_size+1) * sizeof(MPI_Request));
		MPI_Isend(sendbuf, sendcount+rest, sendtype, target, 0, comm, &global_reqs[global_reqs_counter]);
		int i;
		for (i=original_comm_size; i<size; i++) {
			MPI_Isend(&sendbuf[sendcount*(i-original_comm_size+1)+rest], sendcount, sendtype, i, 0, comm, &global_reqs[global_reqs_counter+i-original_comm_size+1]);
		}
		global_reqs_counter += i-original_comm_size+1;
	}
	else if (rank == target) {
		//MPI_Request reqs[size-original_comm_size+1];
		global_reqs = (MPI_Request *) realloc(global_reqs,
				(global_reqs_counter+size-original_comm_size+1) * sizeof(MPI_Request));
		MPI_Irecv(recvbuf, recvcount+rest, recvtype, source, 0, comm, &global_reqs[global_reqs_counter]);
		int i;
		for (i=original_comm_size; i<size; i++) {
			MPI_Irecv(&recvbuf[recvcount*(i-original_comm_size+1)+rest], recvcount, recvtype, i, 0, comm, &global_reqs[global_reqs_counter+i-original_comm_size+1]);
		}
		global_reqs_counter += i-original_comm_size+1;
	}
	// All spawned processes simply forward chunks of data
	else if (rank >= original_comm_size) {
		MPI_Recv(interbuf, recvcount, recvtype, source, 0, comm,
				MPI_STATUS_IGNORE);
		MPI_Send(interbuf, sendcount, sendtype, target, 0, comm);
	}
	//Don't forget to complete communication calls for blocking cases at source and target
	if ((rank == source) || (rank == target)) {
		if (MPIB_p2p_blocking) {
			MPI_Waitall(global_reqs_counter, global_reqs, MPI_STATUSES_IGNORE);
			free(global_reqs);
			global_reqs = NULL;
			global_reqs_counter = 0;
			
		}
	}

	return MPI_SUCCESS;
}

/*!
 * Allocate temporary buffer, call modified point-to-point and
 * release the buffer
 */
int MPIB_Send_sg(void* buf, int count, MPI_Datatype datatype, int dest,
		int tag, MPI_Comm comm, int MPIB_p2p_blocking) {
	int size, rank;
	MPI_Comm_size(intracomm, &size);
	MPI_Comm_rank(intracomm, &rank);
	int num_pieces = size-original_comm_size+2-1;
	int diff = count-(count/num_pieces)*num_pieces;
	MPIB_scatter_gather_based_p2p(buf, (count / num_pieces), diff, datatype,
			NULL, (count/num_pieces), datatype, intracomm,
			rank, dest, NULL, MPIB_p2p_blocking);
	return 0;
}

/*!
 * Allocate temporary buffer, call modified point-to-point and
 * release the buffer
 */
int MPIB_Recv_sg(void* buf, int count, MPI_Datatype datatype, int source,
		int tag, MPI_Comm comm, MPI_Status* status, int MPIB_p2p_blocking) {
	int size, rank;
	MPI_Comm_size(intracomm, &size);
	MPI_Comm_rank(intracomm, &rank);
	int num_pieces = size-original_comm_size+2-1;
	int diff = count-(count/num_pieces)*num_pieces;
	MPIB_scatter_gather_based_p2p(NULL, (count/num_pieces), diff, datatype,
			buf, (count/num_pieces), datatype, intracomm,
			source, rank, NULL, MPIB_p2p_blocking);
	return 0;
}
