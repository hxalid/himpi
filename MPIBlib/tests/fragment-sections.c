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
        int provided;
	double totalTime=0.0;
        int required = MPI_THREAD_MULTIPLE;
        char * empty = NULL;
	int THREAD_NUM;
	THREAD_NUM = 2; 
	omp_set_num_threads(THREAD_NUM);

        MPI_Init_thread(&argc, &argv, required, &provided);
        comm = MPI_COMM_WORLD;
        if (provided != required) {
                fprintf(stderr, "no required provided\n");
                MPI_Abort(comm, -1);
        }
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	array = (char *) malloc(SIZE * sizeof(char));
	if (rank == 0) {
		for (k=0;k < MAX_STEPS; k++) {
			for (j=1; j < MAX_ITER; j++) {
				start[j] = MPI_Wtime();
				#pragma omp parallel 
				{
					#pragma omp sections
					{
						#pragma omp section
						MPI_Send(array, (SIZE*(k+1))/(2*MAX_STEPS), MPI_CHAR, 1, 0, comm);
						#pragma omp section
						MPI_Send(array+SIZE/2, (SIZE*(k+1))/(2*MAX_STEPS), MPI_CHAR, 1, 0, comm);
					}
				}

				MPI_Recv(empty, 0, MPI_CHAR, 1, 0, comm, MPI_STATUS_IGNORE);
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
				#pragma omp parallel
				{
					#pragma omp sections
					{
						#pragma omp section
						MPI_Recv(array, (SIZE*(k+1))/(2*MAX_STEPS), MPI_CHAR, 0, 0, comm, MPI_STATUS_IGNORE);
						#pragma omp section
						MPI_Recv(array+SIZE/2, (SIZE*(k+1))/(2*MAX_STEPS), MPI_CHAR, 0, 0, comm, MPI_STATUS_IGNORE);
					}
				}
				MPI_Send(empty, 0, MPI_CHAR, 0, 0, comm);
			}
		}
	}

	if (rank == 0) {
	}
        MPI_Finalize();
        return 0;
}

