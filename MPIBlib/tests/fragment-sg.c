#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#define MAX_ITER 16
#define MAX_STEPS 16

int main(int argc, char **argv) {

	int SIZE = 1000000;
        char * s_array = NULL;
	char * r_array = NULL;
        double start[MAX_ITER], finish[MAX_ITER];
        MPI_Comm comm;
	MPI_Comm comm_no_last;
	MPI_Comm comm_no_first;
        int rank, size;
	int j,k;
	double totalTime=0.0;
        char * empty = NULL;
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

	s_array = (char *) malloc(SIZE * sizeof(char));

	int color;
	if (rank == (size-1))
		color = MPI_UNDEFINED;
	else
		color = 1;
		
	MPI_Comm_split(comm, color, rank, &comm_no_last);

	if (rank == 0)
		color = MPI_UNDEFINED;
	else
		color = 1;

	MPI_Comm_split(comm, color, rank, &comm_no_first);

	if ((rank < size -1) && (rank > 0))
		intermed_array = (char *) malloc((SIZE/size) * sizeof(char));

	if (rank == 0)
		s_array = (char *) malloc(SIZE * sizeof(char));

	if (rank == (size-1))
		r_array = (char *) malloc(SIZE * sizeof(char));

	for (k=0;k < MAX_STEPS; k++) {
		for (j=1; j < MAX_ITER; j++) {
			if (rank == 0) 
				start[j] = MPI_Wtime();
			
			//everyone gets a piece at r_array
			if (comm_no_last != MPI_COMM_NULL) {
				int new_first;
				MPI_Comm_rank(comm_no_last, &new_first);
				MPI_Scatter(s_array, SIZE/size, MPI_CHAR, intermed_array, SIZE/size, MPI_CHAR, 0, comm_no_last);
			}
			if (comm_no_first != MPI_COMM_NULL) {
				int new_last;
				MPI_Comm_rank(comm_no_first, &new_last);
				MPI_Gather(intermed_array, SIZE/size, MPI_CHAR, r_array, SIZE/size, MPI_CHAR, size-1, comm_no_first);
			}
			if (rank == size-1)
				MPI_Send(empty, 0, MPI_CHAR, 0, 0, comm);
			
			if (rank == 0) {
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

        MPI_Finalize();
        return 0;
}

