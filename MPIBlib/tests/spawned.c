#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "p2p/mpib_p2p.h"

#define MAX_ITER 16
#define MAX_STEPS 32

int main(int argc, char **argv) {

	int j,k;
	double totalTime=0.0;
	int SIZE = 10000000;
        char * s_array = NULL;
	char * r_array = NULL;
        MPI_Comm intercomm, intracomm;
        int rank, size;
	int real_length;
	char *intermed_array = NULL;
	char name[128];
        MPI_Init(&argc, &argv);
	
	MPI_Comm_get_parent(&intercomm);


	if (intercomm != MPI_COMM_NULL) {
		MPI_Intercomm_merge(intercomm, 1, &intracomm);
		MPI_Comm_size(intracomm, &size);
		MPI_Comm_rank(intracomm, &rank);
		SIZE = (SIZE/size) * size;
		MPI_Get_processor_name(name, &real_length);
		printf("I am a spawned process, rank %d of %d, name is %s\n", rank, size, name);


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
			}
		}
		MPI_Comm_free(&intracomm);
		MPI_Comm_free(&intercomm);
	}

	
        MPI_Finalize();
        return 0;
}

