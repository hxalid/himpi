#include "mpib_p2p.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

//This variable holds the original size of MPI_COMM_WORLD before spawning additional processes
int original_comm_size = 0;

// 0 means standard p2p, 1 means modified p2p
int MPIB_p2p_type = 0;

//This fragment size can be used optionally to fragment any messages inside this benchmark
int MPIB_max_fragment_size = INT_MAX;

// Default number of processes to spawn at each node for modified p2p
int PROC_NUM_SPAWNED = 2;

// Used only to create intracomm
MPI_Comm intercomm = MPI_COMM_NULL;

// Global communicator containing both the initial processes
// as well as all spawned processes
MPI_Comm intracomm = MPI_COMM_NULL;

// Used for control messages between any of the initial processes and any of the spawned processes
// This includes information about the type of communication or a termination signal
MPI_Comm control_intracomm = MPI_COMM_NULL;

// This array  holds the names of the machine hosting each rank:
// i.e. all_names[i*128] = hostname of machine where process i runs
char * all_names;

// Used for non-blocking point-to-point and waitall routines.
MPI_Request* global_reqs = NULL;

// Counter for global_reqs
int global_reqs_counter = 0;

int p2p_init_max_fragment_size() {
	char * max_fragment_size_string = getenv("FRAGMENT_SIZE");
	if (NULL != max_fragment_size_string) {
		MPIB_max_fragment_size = atoi(max_fragment_size_string);
		fprintf(stdout, "#fragment size set to %d bytes\n", MPIB_max_fragment_size);
	}
	return 0;
}

static int cmpr(const void *a, const void *b) { 
 return strcmp(*(char **)a, *(char **)b);
}

int p2p_init(MPI_Comm comm, int proc_num_spawned) {
	//set global variable to always use modified p2p
	MPIB_p2p_type = 1;
	PROC_NUM_SPAWNED = proc_num_spawned;
	char spawned_prog[64];
	strcpy(spawned_prog, "./p2p_forward");
	MPI_Comm_get_parent(&intercomm);
	int merge_flag;
	if (intercomm == MPI_COMM_NULL) {
		merge_flag = 0;
		int size;
		MPI_Comm_size(comm, &size);
		original_comm_size = size;
		int err_codes[size * PROC_NUM_SPAWNED];
		MPI_Comm_spawn(spawned_prog, MPI_ARGV_NULL, size * PROC_NUM_SPAWNED,
				MPI_INFO_NULL/*info*/, 0, comm, &intercomm, err_codes);
	} else {
		merge_flag = 1;
	}
	MPI_Intercomm_merge(intercomm, merge_flag, &intracomm);
	int intercomm_size;
	MPI_Comm_size(intercomm, &intercomm_size);
	MPI_Comm_dup(intracomm, &control_intracomm);
	MPI_Bcast(&PROC_NUM_SPAWNED, 1, MPI_INT, 0/*root*/, intracomm);

	int real_length;
	char name[128];
	MPI_Get_processor_name(name, &real_length);
	int intracomm_size, rank;
	MPI_Comm_size(intracomm, &intracomm_size);
	MPI_Comm_rank(intracomm, &rank);
	
	//let spawned processes know about the original comm size
	MPI_Bcast(&original_comm_size, 1, MPI_INT, 0, intracomm);
	all_names = (char *) malloc(intracomm_size * 128 * sizeof(char));
	MPI_Allgather(name, 128, MPI_CHAR, all_names, 128, MPI_CHAR, intracomm);

	int new_rank;	
	int belong_to_source;

	//initial procs don't participate in re-mapping
	if (rank < original_comm_size) {
		belong_to_source = MPI_UNDEFINED;
	}
	//spawned procs participate in re-mapping:
	else {
		//here, we should be careful: we might have experiments
		//on one single machine); in this case, put all nodes
		//on the same node (2 even groups)
		if (strncmp(all_names,all_names+128,128) == 0) {
			belong_to_source = rank % 2;
		}
		//if nodes in different domains, figure a 50/50 mapping
		//(now silly check if there is a DOT in sender's name only_
		else if (strchr(all_names,'.') != NULL) {
			//get pointer to domain name (after '.')
			char * domain = strchr(all_names,'.')+1;
			//domain name of sender is in my name!
			if (strstr(name, domain) != NULL)
				belong_to_source = 1;
			//domain name is not in my name!
			else
				belong_to_source = 0;
		}
		//assume that it's two different nodes on same domain, so just map
		//evenly the helper processes on the other nodes
		else {
			belong_to_source = rank % 2;
		}
	}
	//comm2 is just a helper for generating sub communicators to use
	//for rank mapping
	MPI_Comm comm2;
	MPI_Comm_split(intracomm, belong_to_source,rank, &comm2);
	if (comm2 == MPI_COMM_NULL) printf("DEBUG: comm null for rank %d\n", rank);
	int key;
	if (rank < original_comm_size) {
		new_rank = rank;
		key = 1;
	}
	else {
		MPI_Comm_rank(comm2, &new_rank);
		new_rank += original_comm_size;
		key = 1;
		if (!belong_to_source) {new_rank += (original_comm_size*PROC_NUM_SPAWNED )/2; }
	}
	MPI_Comm_split(control_intracomm, key, new_rank, &intracomm);	
	return 0;
}

int signal_spawned_procs(int rank_list_size, int flag, int source, int target, int sendcount, int rest);

int p2p_finalize() {
	int size, initial_size, rank;
	if (intracomm != MPI_COMM_NULL) {
		MPI_Comm_rank(intracomm, &rank);
		MPI_Comm_size(intracomm, &size);
		if (rank == 0) {
			MPI_Comm comm = MPI_COMM_WORLD;
			MPI_Comm_size(comm, &initial_size);
			signal_spawned_procs(size, 2/*terminate flag*/, 0, 0, 0, 0);
		}
		MPI_Comm_free(&intracomm);
		MPI_Comm_free(&control_intracomm);
		MPI_Comm_free(&intercomm);
	}
	free(all_names);
	return 0;
}

int MPIB_Send_sg(void *buf, int count, MPI_Datatype datatype, int dest,
		int tag, MPI_Comm comm, int MPIB_p2p_blocking);

int MPIB_Send(void *buf, int count, MPI_Datatype datatype, int dest,
		int tag, MPI_Comm comm) {
	if (MPIB_p2p_type == 0) {
		if (count <= MPIB_max_fragment_size)
			MPI_Send(buf, count, datatype, dest, tag, comm);
			//MPI_Rsend(buf, count, datatype, dest, tag, comm);
		//fragmentation of message
		else {
			int fragment;
			for (fragment = MPIB_max_fragment_size; fragment < count; fragment += MPIB_max_fragment_size) {
				MPI_Send(buf, MPIB_max_fragment_size, datatype, dest, tag, comm);
			}
			fragment -= MPIB_max_fragment_size;
			MPI_Send(buf, count-fragment, datatype, dest, tag, comm);

			for (fragment = MPIB_max_fragment_size; fragment < count; fragment += MPIB_max_fragment_size) {
				MPI_Recv(buf, MPIB_max_fragment_size, datatype, dest, tag, comm, MPI_STATUS_IGNORE);
			}
			fragment -= MPIB_max_fragment_size;
			MPI_Recv(buf, count-fragment, datatype, dest, tag, comm, MPI_STATUS_IGNORE);
		}
		return 0;
	}
	else if (MPIB_p2p_type == 1) {
		return MPIB_Send_sg(buf, count, datatype, dest, tag, comm, 1);
	}
	fprintf(stderr, "Wrong MPIB_p2p_type = %d\n", MPIB_p2p_type);
	return -1;

}

int MPIB_Isend(void *buf, int count, MPI_Datatype datatype, int dest,
		int tag, MPI_Comm comm, MPI_Request *request) {
	if (MPIB_p2p_type == 0)
		return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
	else if (MPIB_p2p_type == 1) {
		*request = MPI_REQUEST_NULL;
		return MPIB_Send_sg(buf, count, datatype, dest, tag, comm, 0);
	}
	fprintf(stderr, "Wrong MPIB_p2p_type = %d\n", MPIB_p2p_type);
	return -1;
}

int MPIB_Recv_sg(void *buf, int count, MPI_Datatype datatype, int source,
		int tag, MPI_Comm comm, MPI_Status *status, int MPIB_p2p_blocking);

int MPIB_Recv(void *buf, int count, MPI_Datatype datatype, int source,
		int tag, MPI_Comm comm, MPI_Status *status) {
	if (MPIB_p2p_type == 0) {
		if (count <= MPIB_max_fragment_size) 
			MPI_Recv(buf, count, datatype, source, tag, comm, status);
		else { 
			int fragment;
			for (fragment = MPIB_max_fragment_size; fragment < count; fragment += MPIB_max_fragment_size) {
				MPI_Recv(buf, MPIB_max_fragment_size, datatype, source, tag, comm, MPI_STATUS_IGNORE);
			}
			fragment -= MPIB_max_fragment_size;
			MPI_Recv(buf, count-fragment, datatype, source, tag, comm, MPI_STATUS_IGNORE);
			
			for (fragment = MPIB_max_fragment_size; fragment < count; fragment += MPIB_max_fragment_size) {
				MPI_Send(buf, MPIB_max_fragment_size, datatype, source, tag, comm);
			}
			fragment -= MPIB_max_fragment_size;
			MPI_Send(buf, count-fragment, datatype, source, tag, comm);
		}
		return 0;
	}
	else if (MPIB_p2p_type == 1)
		return MPIB_Recv_sg(buf, count, datatype, source, tag, comm, status, 1);
	fprintf(stderr, "Wrong MPIB_p2p_type = %d\n", MPIB_p2p_type);
	return -1;
}

int MPIB_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
		int tag, MPI_Comm comm, MPI_Request *request) {
	if (MPIB_p2p_type == 0)
		return MPI_Irecv(buf, count, datatype, source, tag, comm, request);
	else if (MPIB_p2p_type == 1) {
		*request = MPI_REQUEST_NULL;
		return MPIB_Recv_sg(buf, count, datatype, source, tag, comm,
				MPI_STATUS_IGNORE, 0);
	}
	fprintf(stderr, "Wrong MPIB_p2p_type = %d\n", MPIB_p2p_type);
	return -1;
}

int MPIB_Waitall(int count, MPI_Request *array_of_requests,
		MPI_Status *array_of_statuses) {
	if (MPIB_p2p_type == 0) {
		return MPI_Waitall(count, array_of_requests, array_of_statuses);
	}
	//if we are doing our own p2p communication, then we ignore all the input arguments.
	//Instead, we use our internal request list
	else if (MPIB_p2p_type == 1) {
		int retValue = MPI_Waitall(global_reqs_counter, global_reqs, MPI_STATUSES_IGNORE);
		free(global_reqs);
		global_reqs = NULL;
		global_reqs_counter = 0;
		return retValue;
	}
	fprintf(stderr, "Wrong MPIB_p2p_type = %d\n", MPIB_p2p_type);
	return -1;
}
