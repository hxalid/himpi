#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "p2p/mpib_p2p.h"

#define MAX_ITER 16
#define MAX_STEPS 32
#define PROC_NUM_SPAWNED 4

int main(int argc, char **argv) {

	int i;
        char * s_array = NULL;
	char * r_array = NULL;
	int debug_counter = 0;
        MPI_Comm comm;
        int rank, size;
	int j,k;
	double totalTime=0.0;
	double minTime;
	char name[128];
	int real_length;

	MPI_Init(&argc, &argv);
        comm = MPI_COMM_WORLD;
        p2p_init(comm, PROC_NUM_SPAWNED);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);


	//reassigning size/rank for new intracommunicator
	MPI_Comm_size(intracomm, &size);
	MPI_Comm_rank(intracomm, &rank);
	MPI_Get_processor_name(name, &real_length);
	printf("I am  rank %d of %d, name is %s\n", rank, size, name);


	//rank 0 initializes the data
	if (rank == 0) {
		s_array = (char *) malloc(SIZE * sizeof(char));
		for (i=0;i<SIZE;i++) s_array[i] = 'g';
	}

	if (rank == 1)
		r_array = (char *) malloc(SIZE * sizeof(char));



	for (k=0;k < MAX_STEPS; k++) {
		totalTime = 0.0;
		minTime = 10000.;
		int current_size = (SIZE * (k+1))/MAX_STEPS;
		for (j=0; j < MAX_ITER; j++) {
			if ((rank == 1) && (j == 0)) {
				for (i=0; i<SIZE; i++) 
					r_array[i] = 'X'; 
			}
			MPI_Barrier(comm);
			double start = MPI_Wtime();
			if (rank == 0)
				MPIB_Send_sg(s_array, current_size, MPI_CHAR, 1, 0, comm);
			if (rank == 1)
				MPIB_Recv_sg(r_array, current_size, MPI_CHAR, 0, 0, comm, MPI_STATUS_IGNORE);
			double end = MPI_Wtime();
			if ((rank ==1) && (j == 0)) {
				debug_counter=0;
				for (i=0; i<SIZE; i++) {
					if (r_array[i] == 'g') {
						debug_counter++;
					}
				}
			}
			double duration = end - start;
			double max_duration;
			MPI_Reduce(&duration, &max_duration, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
			if (rank == 0) {
				totalTime+=max_duration;
				if (max_duration < minTime)
					minTime = max_duration;
			}
		}
		/*
		if (rank == 0) {
			printf("%d\t%lf\n", current_size,  minTime);
		}
		*/
		MPI_Bcast(&debug_counter, 1, MPI_INT, 1, comm);
		if (rank == 0) {
			printf("%d\t%lf\n", debug_counter, minTime);
		}
	}

	p2p_finalize();
        MPI_Finalize();
        return 0;
}

