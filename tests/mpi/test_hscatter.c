#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include <getopt.h>  

#define MAX_INT_SIZE 50
#define MIN_PROCS 4

typedef struct {
	int rec_out_group_comm;
	int rec_in_group_comm;
	int rec_comm_world;
	int err_flag;
	double in_time;
	double out_time;
	double group_time;
	double no_group_time;
} t_mpi_coll_response;

void print_array(int my_rank, int size, double* array) {
	int i = 0;
	for (i = 0; i < size; i++) {
		printf("[%d]:array[%d]=%f ", my_rank, i, array[i]);
	}
	printf("\n");
}

int *create_rand_elms(int num_elements) {
	int *rand_nums = (int *) malloc(sizeof(int) * num_elements);
	assert(rand_nums != NULL);
	int i;
	for (i = 0; i < num_elements; i++) {
		rand_nums[i] = (rand() / (int) MAX_INT_SIZE);
	}
	return rand_nums;
}

// Computes the average of an array of numbers
double compute_avg(double *array, int num_elements) {
	double sum = 0.0;
	int i;
	for (i = 0; i < num_elements; i++) {
		sum += array[i];
	}
	return sum / num_elements;
}

int MPI_HGather(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm, int num_groups) {
	if (num_groups == 1) {
		return MPI_Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt,
				recvtype, root, comm);
	} else {
		int num_procs;
		int rank;
		MPI_Comm in_group_comm, out_group_comm;

		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &num_procs);

		int pg = num_procs / num_groups;
		int my_group = rank / pg;

		if (num_procs % num_groups != 0 || num_procs <= num_groups) {
			fprintf(stderr, "MPI_HGather: Wrong number of groups\n");
			return -1;
		}

		int data_type_size;
		MPI_Type_size(recvtype, &data_type_size);

		void* g_rcvbuf = malloc(data_type_size * pg * recvcnt);

		MPI_Comm_split(comm, (rank - my_group * pg == 0) ? 0 : MPI_UNDEFINED,
				rank, &out_group_comm);
		MPI_Comm_split(comm, my_group, rank, &in_group_comm);

		/* Gather inside groups */
		if (in_group_comm != MPI_COMM_NULL)
			MPI_Gather(sendbuf, recvcnt, sendtype, g_rcvbuf, recvcnt, sendtype,
					root, in_group_comm);

		/* Gather between groups */
		if (out_group_comm != MPI_COMM_NULL)
			MPI_Gather(g_rcvbuf, pg * recvcnt, sendtype, recvbuf, pg * recvcnt,
					sendtype, root, out_group_comm);

		if (out_group_comm != MPI_COMM_NULL)
			MPI_Comm_free(&out_group_comm);
		if (in_group_comm != MPI_COMM_NULL)
			MPI_Comm_free(&in_group_comm);

		free(g_rcvbuf);

		return 0;
	}
}

/*
 make the send and receive buffer types general
 */
int MPI_HScatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
		MPI_Comm comm, int num_groups) {
	if (num_groups == 1) {
		return MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt,
				recvtype, root, comm);
	} else {
		int num_procs;
		int rank;
		MPI_Comm in_group_comm, out_group_comm;

		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &num_procs);

		if (num_procs % num_groups != 0 || num_procs <= num_groups) {
			fprintf(stderr, "HScatter: Wrong number of groups\n");
			return -1;
		}

		int pg = num_procs / num_groups;
		int my_group = rank / pg;

		int data_type_size;
		MPI_Type_size(recvtype, &data_type_size);

		void* g_rcvbuf = malloc(data_type_size * pg * recvcnt);

		MPI_Comm_split(comm, (rank - my_group * pg == 0) ? 0 : MPI_UNDEFINED,
				rank, &out_group_comm);
		MPI_Comm_split(comm, my_group, rank, &in_group_comm);

		/* Scatter between groups */
		if (out_group_comm != MPI_COMM_NULL) {
			MPI_Scatter(sendbuf, pg * recvcnt, sendtype, g_rcvbuf, pg * recvcnt,
					recvtype, root, out_group_comm);
			MPI_Comm_free(&out_group_comm);
		}

		/* Scatter inside groups */
		if (in_group_comm != MPI_COMM_NULL) {
			MPI_Scatter(g_rcvbuf, recvcnt, sendtype, recvbuf, recvcnt, recvtype,
					root, in_group_comm);
			MPI_Comm_free(&in_group_comm);
		}

		free(g_rcvbuf);

		return 0;
	}
}

int main(int argc, char** argv) {
	/*
	 *  These values are used if not specified by the user 
	 */
	int i = 0;
	int rec = 1;
	int root = 0;
	int reps = 30;
	int debug = 0;
	int use_any_root = 0;
	int dont_use_comm_world = 0;
//    int pg = 1;
	int groups = 1;
//    int my_group = 0;
	int rank;
	int num_proc;
	int num_proc_new;
	int msg_factor = 2;
	int msg_size;
	int rec_in_group_comm = -1;
	int rec_out_group_comm = -1;
	int rec_comm_world = -1;
	int p_name_len;
	char p_name[MPI_MAX_PROCESSOR_NAME];

	int* snd_array; /* storage for message  */
	int* rcv_array;

	double start_time, elapsed_time, max_time; /* time spent for collectives */

	srand(time(NULL));

	/* Start up MPI */
	MPI_Init(&argc, &argv);

	/* Find out processor name */
	//   MPI_Get_processor_name(p_name, &p_name_len);
	/* Find out process rank  */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* Find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	int g_max = num_proc;

	int opt;
	while ((opt = getopt(argc, argv, "m:r:s:i:a:g:f:d:wuh")) != -1) {
		switch (opt) {
		case 'm':
			msg_size = atoi(optarg);
			break;
		case 'r':
			rec = atoi(optarg);
			break;
		case 's':
			root = atoi(optarg);
			break;
		case 'i':
			reps = atoi(optarg);
			break;
		case 'g':
			g_max = atoi(optarg);
			break;
		case 'f':
			msg_factor = atoi(optarg);
			break;
		case 'd':
			debug = atoi(optarg);
			break;
		case 'u':
			use_any_root = 1;
			break;
		case 'w':
			dont_use_comm_world = 1;
			break;
		case 'h':
			if (rank == 0)
				fprintf(stdout, "%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
						"-m message size", "-f message_factor",
						"-r number of hierarchies", "-s scatter root",
						"-i number of iterations",
						"-g maximum number of groups",
						"-d set to 0 disable host names. Default is 1",
						"-w if used don't create sub-communicators", "-h help");
			MPI_Finalize();
			return 0;
		default:
			fprintf(stderr, "Unknown option %s\n", optarg);
			break;
		}
	}

	if ((use_any_root && (root < 0 || root > num_proc))
			|| (!use_any_root && root)) {
		if (!rank)
			fprintf(stdout, "Wrong root\n");
		MPI_Finalize();
		return -1;
	}

	MPI_Comm comm_measure;
	if (dont_use_comm_world) {
		MPI_Comm_split(MPI_COMM_WORLD, (rank > 1) ? 0 : MPI_UNDEFINED, rank,
				&comm_measure);
	} else {
		comm_measure = MPI_COMM_WORLD;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (comm_measure != MPI_COMM_NULL) {

		/* Find out processor name */
		MPI_Get_processor_name(p_name, &p_name_len);

		/* Find out process rank  */
		MPI_Comm_rank(comm_measure, &rank);

		/* Find out number of processes */
		MPI_Comm_size(comm_measure, &num_proc_new);

		g_max = num_proc_new;

		if (debug > 0) {
			fprintf(stdout, "%d: %s\n", rank, p_name);
			MPI_Barrier(comm_measure);
			MPI_Barrier(comm_measure);
		}

		if (rank == root) {
			fprintf(stdout,
					"Inputs: p=%d, p_new=%d, root=%d, m=%d, f=%d, rec=%d, reps=%d, d=%d, dont_use_comm_world=%d, rec_comm_world=%d\n",
					num_proc, num_proc_new, root, msg_size, msg_factor, rec,
					reps, debug, dont_use_comm_world, rec_comm_world);

			fprintf(stdout,
					"root  p   g    msg(KB)     time     rec   reps   rec_in  rec_out  rec_world\n");
		}

		/*
		 * Start benchmarking hscatter
		 */
		if (rec != 0) {
			groups = 1;
			int hscatter_res;
			while (groups < g_max) {
				if (num_proc_new % groups == 0) {

					if (rank == root) {
						//msg_size is number of elements per process
						snd_array = create_rand_elms(msg_size); // on the root we have very large array
					}
					rcv_array = (int*) malloc(sizeof(int) * msg_size / num_proc_new);
					assert(rcv_array != NULL);

					MPI_Barrier(comm_measure);

					elapsed_time = 0.;
					for (i = 0; i < reps; i++) {
						MPI_Barrier(comm_measure);

						start_time = MPI_Wtime();

						hscatter_res = MPI_HScatter(snd_array, msg_size,
								MPI_INT, rcv_array, msg_size, MPI_INT, root,
								comm_measure, groups);

						if (hscatter_res == -1) {
							MPI_Abort(comm_measure, 1986);
							MPI_Finalize();
						}

						elapsed_time += (MPI_Wtime() - start_time);

						MPI_Barrier(comm_measure);
					}

					elapsed_time /= reps;
					MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX,
							root, comm_measure);

					if (rank == root) {
						fprintf(stdout, "%d %d %d %f %f %d %d %d %d %d\n", root,
								num_proc_new, groups,
								msg_size * sizeof(int) / 1000., max_time, rec,
								reps, rec_in_group_comm, rec_out_group_comm,
								rec_comm_world);
					}

					MPI_Barrier(comm_measure);

					if (rank == root)
						free(snd_array);
					free(rcv_array);

				}   // if (num_proc_new % groups == 0)

				groups++;
			}  // while

		}   // if (rec != 0)

		// Measuring comm split time
		/*
		 if (debug == 2) {
		 MPI_Comm in_group_comm, out_group_comm;
		 double max_time_out, max_time_in;

		 if (num_proc_new > MIN_PROCS) {
		 if (rank == root)
		 fprintf(stdout, "SPLIT. root p g max_time_in max_time_out rec reps\n");

		 groups = 2;

		 while (groups < g_max) {
		 if ( num_proc_new % groups == 0 ) {
		 pg = num_proc_new / groups;
		 my_group = rank / pg;

		 elapsed_time = 0.;
		 for (i = 0; i < reps; i++) {
		 MPI_Barrier(comm_measure);
		 start_time = MPI_Wtime();
		 MPI_Comm_split(comm_measure, (rank - my_group * pg == 0) ? 0 : MPI_UNDEFINED, rank, &out_group_comm);
		 elapsed_time += (MPI_Wtime() - start_time);
		 }
		 elapsed_time /= reps;
		 MPI_Barrier(comm_measure);
		 MPI_Reduce(&elapsed_time, &max_time_out, 1, MPI_DOUBLE, MPI_MAX, root, comm_measure);

		 elapsed_time = 0.;
		 for (i = 0; i < reps; i++) {
		 MPI_Barrier(comm_measure);
		 start_time = MPI_Wtime();
		 MPI_Comm_split(comm_measure, my_group, rank, &in_group_comm);
		 elapsed_time += (MPI_Wtime() - start_time);
		 }
		 elapsed_time /= reps;
		 MPI_Barrier(comm_measure);
		 MPI_Reduce(&elapsed_time, &max_time_in, 1, MPI_DOUBLE, MPI_MAX, root, comm_measure);

		 if (rank == root) {
		 fprintf(stdout, "%d %d %d %f %f %d %d\n",
		 root,
		 num_proc_new,
		 groups,
		 max_time_in,
		 max_time_out,
		 rec,
		 reps
		 );
		 }

		 if (out_group_comm!=MPI_COMM_NULL)
		 MPI_Comm_free(&out_group_comm);
		 if (in_group_comm!=MPI_COMM_NULL)
		 MPI_Comm_free(&in_group_comm);
		 }
		 groups++;
		 }
		 }
		 }

		 */

	}  // if (comm_measure != MPI_COMM_NULL )

	MPI_Barrier(MPI_COMM_WORLD);
	/* Shut down MPI */
	MPI_Finalize();

	return EXIT_SUCCESS;
}
