#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>

#define MAX_ITER 16
#define MAX_STEPS 32

int MPIB_Gather_flat(void* sendbuf, int sendcount, MPI_Datatype sendtype, 
	void* recvbuf, int recvcount, MPI_Datatype recvtype,
	int root, MPI_Comm comm)
{
	int rank;
	MPI_Comm_rank(comm, &rank);
	if (rank == root) {
		int size;
		MPI_Comm_size(comm, &size);
		MPI_Aint extent;
		MPI_Type_extent(recvtype, &extent);
		int inc = recvcount * extent;
		int i;
		char* ptr;
		for (i = 0, ptr = (char*)recvbuf; i < size; i++, ptr += inc) {
			if (i == root)
				memcpy(ptr, sendbuf, inc);
			else
				MPI_Recv(ptr, recvcount, recvtype, i, 0, comm, MPI_STATUS_IGNORE);
		}
	} else
		MPI_Send(sendbuf, sendcount, sendtype, root, 0, comm);
	return MPI_SUCCESS;
}

int MPIB_Scatter_flat(void* sendbuf, int sendcount, MPI_Datatype sendtype, 
	void* recvbuf, int recvcount, MPI_Datatype recvtype, 
	int root, MPI_Comm comm)
{
	int rank;
	MPI_Comm_rank(comm, &rank);
	if (rank == root) {
		int size;
		MPI_Comm_size(comm, &size);
		MPI_Aint extent;
		MPI_Type_extent(sendtype, &extent);
		int inc = sendcount * extent;
		int i;
		char* ptr;
		for (i = 0, ptr = (char*)sendbuf; i < size; i++, ptr += inc) {
			if (i == root) {
				MPI_Aint extent;
				MPI_Type_extent(recvtype, &extent);
				memcpy(recvbuf, ptr, recvcount * extent);
			} else
				MPI_Send(ptr, sendcount, sendtype, i, 0, comm);
		}
	} else
		MPI_Recv(recvbuf, recvcount, recvtype, root, 0, comm, MPI_STATUS_IGNORE);
	return MPI_SUCCESS;
}

int main(int argc, char **argv) {

	int i;
	int SIZE = 1000000;
        char * s_array = NULL;
	char * r_array = NULL;
        MPI_Comm comm;
        int rank, size;
	int j,k;
	double totalTime=0.0;
	char name[128];
	char *intermed_array = NULL;
	int real_length;

        MPI_Init(&argc, &argv);
        comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	SIZE = (SIZE/size) * size;
	MPI_Get_processor_name(name, &real_length);
	printf("I am  rank %d, name is %s\n", rank, name);

	intermed_array = (char *) malloc((SIZE/size) * sizeof(char));

	if (rank == 0) {
		s_array = (char *) malloc(SIZE * sizeof(char));
		for (i=0;i<SIZE;i++) s_array[i] = 'g';
	}

	if (rank == (size-1))
		r_array = (char *) malloc(SIZE * sizeof(char));

	for (k=0;k < MAX_STEPS; k++) {
		totalTime = 0.0;
		for (j=0; j < MAX_ITER; j++) {
			MPI_Barrier(comm);
			double start = MPI_Wtime();
			MPIB_Scatter_flat(s_array, (SIZE*(k+1))/(size*MAX_STEPS), MPI_CHAR, intermed_array,(SIZE*(k+1))/(size*MAX_STEPS), MPI_CHAR, 0, comm);
			MPIB_Gather_flat(intermed_array, (SIZE*(k+1))/(size*MAX_STEPS), MPI_CHAR, r_array,SIZE*(k+1)/(size*MAX_STEPS), MPI_CHAR, size-1, comm);
			double end = MPI_Wtime();
			double duration = end - start;
			double max_duration;
			MPI_Reduce(&duration, &max_duration, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
			if (rank == 0) {
					totalTime+=max_duration;
			}
		}
		if (rank == 0) {
			printf("%d\t%lf\n",(SIZE*(k+1))/(MAX_STEPS), totalTime/(double)MAX_ITER);
		}
	}

        MPI_Finalize();
        return 0;
}

