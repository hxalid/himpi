#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#define SIZE 10000000
#define MAX_ITER 16
#define MAX_STEPS 16

int main(int argc, char **argv) {

        char * array;
        double start[MAX_ITER], finish[MAX_ITER];
        MPI_Comm comm;
        int rank, size;
	int j,k;
	double totalTime=0.0;
        char * empty = NULL;
	char name[128];
	int real_length;

        MPI_Init(&argc, &argv);
        comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	MPI_Get_processor_name(name, &real_length);
	printf("I am  rank %d, name is %s\n", rank, name);

	array = (char *) malloc(SIZE * sizeof(char));
	if (rank == 0) {
		for (k=0;k < MAX_STEPS; k++) {
			for (j=1; j < MAX_ITER; j++) {
				start[j] = MPI_Wtime();
				MPI_Request reqs[2];
				//intermediate
				MPI_Isend(array,  (SIZE*(k+1))/(2*MAX_STEPS), MPI_CHAR, 1, 0, comm, &reqs[0]);
				//final receiver
				MPI_Isend(array+SIZE/2,  (SIZE*(k+1))/(2*MAX_STEPS), MPI_CHAR, 2, 0, comm, &reqs[1]);
				MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);
				MPI_Recv(empty, 0, MPI_CHAR, 2, 0, comm, MPI_STATUS_IGNORE);
				finish[j] = MPI_Wtime();
			}
			totalTime=0.0;
			for (j=1; j<MAX_ITER; j++) {
				totalTime+=finish[j]-start[j];
			}
			printf("%d\t%lf\n",(SIZE*(k+1))/(MAX_STEPS), totalTime/(double)MAX_ITER);
			
				
		}
	}

	if (rank == 1) {
		for (k=0;k < MAX_STEPS; k++) {
			for (j=1; j < MAX_ITER; j++) {
				MPI_Recv(array, (SIZE*(k+1))/(2*MAX_STEPS), MPI_CHAR, 0, 0, comm, MPI_STATUS_IGNORE);
				MPI_Send(array, (SIZE*(k+1))/(2*MAX_STEPS), MPI_CHAR, 2, 0, comm);
			}
		}
	}
	if (rank == 3) {
		for (k=0;k < MAX_STEPS; k++) {
			for (j=1; j < MAX_ITER; j++) {
				MPI_Recv(array, (SIZE*(k+1))/(2*MAX_STEPS), MPI_CHAR, 0, 0, comm, MPI_STATUS_IGNORE);
				MPI_Send(array, (SIZE*(k+1))/(2*MAX_STEPS), MPI_CHAR, 2, 0, comm);
			}
		}
	}

	if (rank == 2) {
		for (k=0;k < MAX_STEPS; k++) {
			for (j=1; j < MAX_ITER; j++) {
				MPI_Request reqs[2];
				//from intermediate
				MPI_Irecv(array, (SIZE*(k+1))/(2*MAX_STEPS), MPI_CHAR, 1, 0, comm, &reqs[0]);
				//from original sender
				MPI_Irecv(array+SIZE/2, (SIZE*(k+1))/(2*MAX_STEPS), MPI_CHAR, 0, 0, comm, &reqs[1]);
				MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);
				MPI_Send(empty, 0, MPI_CHAR, 0, 0, comm);
			}
		}
	}

        MPI_Finalize();
        return 0;
}

