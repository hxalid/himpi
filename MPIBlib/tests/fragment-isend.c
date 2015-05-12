#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

#define SIZE 10000000
#define MAX_ITER 16
#define MAX_STEPS 16
#define THREAD_NUM 2

int main(int argc, char **argv) {

        char * array;
        double start[MAX_ITER], finish[MAX_ITER];
        MPI_Comm comm;
        int rank, size;
	int i,j,k;
	double totalTime=0.0;
        char * empty = NULL;

        MPI_Init(&argc, &argv);
        comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	MPI_Request reqs[THREAD_NUM];

	array = (char *) malloc(SIZE * sizeof(char));
	if (rank == 0) {
		for (k=0;k < MAX_STEPS; k++) {
			for (j=1; j < MAX_ITER; j++) {
				start[j] = MPI_Wtime();
				for (i=0; i< THREAD_NUM; i++) {
					MPI_Isend(array+(SIZE*i)/THREAD_NUM, /*SIZE/THREAD_NUM*/ (SIZE*(k+1))/(MAX_STEPS*THREAD_NUM), MPI_CHAR, 1, 0, comm, reqs+i );
				}

				MPI_Recv(empty, 0, MPI_CHAR, 1, 0, comm, MPI_STATUS_IGNORE);
				MPI_Waitall(THREAD_NUM, reqs, MPI_STATUSES_IGNORE);
				finish[j] = MPI_Wtime();
			}
			totalTime=0.0;
			for (j=1; j<MAX_ITER; j++) {
				totalTime+=finish[j]-start[j];
			}
			printf("%d\t%lf\n",(SIZE*(k+1))/(MAX_STEPS*THREAD_NUM), totalTime/(double)MAX_ITER);
			
				
		}
	}

	if (rank == 1) {
		for (k=0;k < MAX_STEPS; k++) {
			for (j=1; j < MAX_ITER; j++) {
				for (i=0; i< THREAD_NUM; i++) {
					MPI_Irecv(array+(SIZE*i)/THREAD_NUM, (SIZE*(k+1))/(MAX_STEPS*THREAD_NUM), MPI_CHAR, 0, 0, comm, reqs+i);
				}
				MPI_Waitall(THREAD_NUM, reqs, MPI_STATUSES_IGNORE);
				MPI_Send(empty, 0, MPI_CHAR, 0, 0, comm);
			}
		}
	}

	if (rank == 0) {
	}
        MPI_Finalize();
        return 0;
}

