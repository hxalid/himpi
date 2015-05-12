#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

#define SIZE 32768 
#define MAX_ITER 5097
#define MAX_STEPS 1

int main(int argc, char **argv) {

        char * array;
        double start[MAX_ITER], finish[MAX_ITER];
        MPI_Comm comm;
        int rank, size;
	int j,k;
	double totalTime=0.0;
	double minTime=1000.;

        MPI_Init(&argc, &argv);
        comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	array = (char *) malloc(SIZE * sizeof(char));
	if (rank == 0) {
		double time[MAX_ITER];
		double finalTime[MAX_ITER];
		for (k=0;k < MAX_STEPS; k++) {
			int current_size = (SIZE*(k+1))/MAX_STEPS;
			start[1] = MPI_Wtime();
			for (j=1; j < MAX_ITER; j++) {
				MPI_Send(array, current_size, MPI_CHAR, 1, 0, comm);
			}
			finish[1] = MPI_Wtime();
			time[1] = finish[1] - start[1];
			MPI_Reduce(&time[1], &finalTime[1], 1, MPI_DOUBLE, MPI_MAX, 0, comm);
			totalTime=0.0;
			minTime = 100000.;
			for (j=1; j<MAX_ITER; j++) {
				minTime= (finalTime[j] < minTime)?finalTime[j]:minTime;
			}
			totalTime+=finalTime[1];
			printf("totalTime  for message of size %d is %lf\n", current_size, totalTime);
			
				
		}
	}

	if (rank == 1) {
		double time[MAX_ITER];
		double finalTime[MAX_ITER];
		for (k=0;k < MAX_STEPS; k++) {
			int current_size = (SIZE*(k+1))/MAX_STEPS;
			start[1] = MPI_Wtime();
			for (j=1; j < MAX_ITER; j++) {
				MPI_Recv(array, current_size, MPI_CHAR, 0, 0, comm, MPI_STATUS_IGNORE);
			}
			finish[1] = MPI_Wtime();
			time[1] = finish[1] - start[1];
			MPI_Reduce(&time[1], &finalTime[1], 1, MPI_DOUBLE, MPI_MAX, 0, comm);
		}
	}

	if (rank == 0) {
	}
        MPI_Finalize();
        return 0;
}

