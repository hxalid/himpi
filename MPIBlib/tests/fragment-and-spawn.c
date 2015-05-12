#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "p2p/mpib_p2p.h"

#define MAX_ITER 16
#define MAX_STEPS 32
#define PROC_NUM_SPAWNED 32


int main(int argc, char **argv) {

	int i;
	int SIZE = 10000000;
        char * s_array = NULL;
	char * r_array = NULL;
        MPI_Comm comm, intercomm, intracomm;
        int rank, size;
	int j,k;
	double totalTime=0.0;
	char name[128];
	char *intermed_array = NULL;
	int real_length;
	int err_codes[PROC_NUM_SPAWNED];

        MPI_Init(&argc, &argv);
        comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	char spawned_prog[64];
	strcpy(spawned_prog, "./spawned");

	MPI_Get_processor_name(name, &real_length);

	//process (size-1) spawns 8 extra processes
	MPI_Bcast(name, 128, MPI_CHAR, 1, comm);
	
	MPI_Info info;
	MPI_Info_create(&info);
	MPI_Info_set(info, "host", name);
	MPI_Comm_spawn(spawned_prog, MPI_ARGV_NULL, PROC_NUM_SPAWNED, info, 1, comm, &intercomm, err_codes);
	MPI_Intercomm_merge(intercomm, 0, &intracomm);


	//reassigning size/rank for new intracommunicator
	MPI_Comm_size(intracomm, &size);
	MPI_Comm_rank(intracomm, &rank);
	MPI_Get_processor_name(name, &real_length);
	printf("I am  rank %d of %d, name is %s\n", rank, size, name);

	SIZE = (SIZE/size) * size;

	//rank 0 initializes the data
	if (rank == 0) {
		s_array = (char *) malloc(SIZE * sizeof(char));
		for (i=0;i<SIZE;i++) s_array[i] = 'g';
	}

	if (rank == 1)
		r_array = (char *) malloc(SIZE * sizeof(char));

	intermed_array = (char *) malloc((SIZE/size) * sizeof(char));


	for (k=0;k < MAX_STEPS; k++) {
		totalTime = 0.0;
		for (j=0; j < MAX_ITER; j++) {
			MPI_Barrier(intracomm);
			double start = MPI_Wtime();
			MPIB_Scatter_flat(s_array, (SIZE*(k+1))/(size*MAX_STEPS), MPI_CHAR, intermed_array,(SIZE*(k+1))/(size*MAX_STEPS), MPI_CHAR, 0, intracomm);
			MPIB_Gather_flat(intermed_array, (SIZE*(k+1))/(size*MAX_STEPS), MPI_CHAR, r_array,SIZE*(k+1)/(size*MAX_STEPS), MPI_CHAR, 1, intracomm);
			double end = MPI_Wtime();
			double duration = end - start;
			double max_duration;
			MPI_Reduce(&duration, &max_duration, 1, MPI_DOUBLE, MPI_MAX, 0, intracomm);
			if (rank == 0) {
					totalTime+=max_duration;
			}
		}
		if (rank == 0) {
			printf("%d\t%lf\n",(SIZE*(k+1))/(MAX_STEPS), totalTime/(double)MAX_ITER);
		}
	}

	MPI_Comm_free(&intracomm);
	MPI_Comm_free(&intercomm);
        MPI_Finalize();
        return 0;
}

