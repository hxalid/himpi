#include <mpi.h>

#define CHUNKS 5096 
#define CHUNK_SIZE 32768
int main(int argc, char **argv) {
int i;
char *data;
MPI_Init(&argc, &argv);
int rank, size;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
	data = (char *) malloc ( CHUNK_SIZE * sizeof(char));
double t1 = MPI_Wtime();
for (i=0; i <  CHUNKS; i ++) 
MPI_Bcast(data, CHUNK_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
double t2 = MPI_Wtime();
t2 = t2 - t1;
double max_time;
MPI_Reduce(&t2, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		printf("TOTAL TIME: %lf\n", max_time);
	}
MPI_Finalize();
return 0;
}
