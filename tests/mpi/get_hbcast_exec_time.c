#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "MPIBlib/benchmarks/mpib.h"
#include "config.h"
#include "himpi/himpi.h"
#include "tools/utils.h"
#include "tools/optimal_groups.h"

int main(int argc, char* argv[]) {
	int i = 0;
	int rec = 1; //number of hierarchy
	int alg = 0; //broadcast algorithm id
	int root = 0;
	int reps = 30;
	int rec_world = -1;
	int rank;
	int num_proc;
	int message_min = 1;
	int message_max = 10485760; //10MB
	int msg_factor = 2;
	int msg_size;
	char* array; /* storage for message  */
	double start_time, end_time;
	double elapsed_time, max_time; /* time spent for broadcast */

	srand(time(NULL));

	/*! Start up MPI */
	HiMPI_Init(&argc, &argv);

	/*! Find out process rank  */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/*! Find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	int g_max = num_proc;

	int opt;

	if (rank == root) {
		fprintf(stdout,
				"Inputs: root=%d, m=%d, M=%d, f=%d, rec=%d, reps=%d, a=%d, d=%d\n",
				root, message_min, message_max, msg_factor, rec, reps, alg,
				debug);
		fprintf(stdout,
				"root  p   g    msg     time     rec   reps   rec_in  rec_out  rec_world \n");
	}

	int count_m = 0;
	for (msg_size = message_min; msg_size < message_max + 1; msg_size *=
			msg_factor) {
		array = create_rand_elms(msg_size);

		elapsed_time = 0.;
		for (i = 0; i < reps; i++) {
			MPI_Barrier(MPI_COMM_WORLD);

			start_time = MPI_Wtime();
			//hpnla_bcast(array, msg_size, MPI_CHAR, root, MPI_COMM_WORLD, alg);
			MPI_Bcast(array, msg_size, MPI_CHAR, root, MPI_COMM_WORLD);
			elapsed_time += (MPI_Wtime() - start_time);
		}
		elapsed_time /= reps; // time in sec

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, root,
				MPI_COMM_WORLD);

		if (rank == root) {
			fprintf(stdout, "MPI_Bcast: %d %d %d %f %f %d %d %d %d %d\n", root,
					num_proc, 1, msg_size * sizeof(char) / 1024., max_time, rec,
					reps, -1, -1, rec_world);
		}

		free(array);
	}


	if (rec != 0) {
		for (msg_size = message_min; msg_size < message_max + 1; msg_size *=
				msg_factor) {
			array = create_rand_elms(msg_size);
			MPI_Barrier(MPI_COMM_WORLD);

			elapsed_time = 0.;
			for (i = 0; i < reps; i++) {
				MPI_Barrier(MPI_COMM_WORLD);
				start_time = MPI_Wtime();

				HiMPI_Bcast(array, msg_size, MPI_CHAR, root, MPI_COMM_WORLD);
				elapsed_time += (MPI_Wtime() - start_time);
			}

			elapsed_time /= reps;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, root,
					MPI_COMM_WORLD);

			if (rank == root) {
				fprintf(stdout, "HMPI_Bcast: %d %d %f %f %d %d %d\n", root,
						num_proc, msg_size * sizeof(char) / 1024., max_time,
						rec, reps, rec_world);
			}

			free(array);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/*! Shut down MPI */
	HiMPI_Finalize();
} /* main */
