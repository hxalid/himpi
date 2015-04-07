#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <mpi.h>

#include "config.h"

//#include "hbcast_old.h"
#include "../hmpi/hmpi.h"
#include "../tools/utils.h"

#ifdef HAVE_MPIX_H
#include <mpix.h>
#endif


int main(int argc, char* argv[]) {
    int i = 0;
    int rec = 1;
    int alg = 0;
    int root = 0;
    int reps = 30;
    int debug = 0;
    int use_any_root = 0;
    int rec_world = -1;
    int pg = 1;
    int groups = 1;
    int my_group = 0;
    int rank;
    int num_proc;
    int message_min = 1;
    int message_max = 10485760; //10MB
    int msg_factor = 2;
    int msg_size;
    char* array; /* storage for message  */
    double start_time, end_time;
    double elapsed_time, max_time; /* time spent for broadcast */
    int bcast_linear_segment_size = 1024; //TODO


    int p_name_len;
    char p_name[MPI_MAX_PROCESSOR_NAME];

    srand(time(NULL));

    /* Start up MPI */
    MPI_Init(&argc, &argv);

    /* Find out processor name */
    MPI_Get_processor_name(p_name, &p_name_len);

    /* Find out process rank  */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    int g_max = num_proc;

    int opt;
    while ((opt = getopt(argc, argv, "m:M:r:s:i:a:g:f:c:d:uh")) != -1) {
        switch (opt) {
            case 'm':
                message_min = atoi(optarg);
                break;
            case 'M':
                message_max = atoi(optarg);
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
            case 'a':
                alg = atoi(optarg);
                break;
            case 'g':
                g_max = atoi(optarg);
                break;
            case 'f':
                msg_factor = atoi(optarg);
                break;
            case 'c':
                bcast_linear_segment_size = atoi(optarg);
                break;
            case 'd':
                debug = atoi(optarg);
                break;
            case 'u':
                use_any_root = 1;
                break;
            case 'h':
                if (rank == 0)
                    fprintf(stdout,
                        "%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
                        "-m minimum message size",
                        "-M maximum message size",
                        "-f message_factor",
                        "-r number of hierarchies",
                        "-s bcast root",
                        "-i number of iterations",
                        "-a broadcast algorithm to use",
                        "-g maximum number of groups",
                        "-d set to 0 disable host names. Default is 1",
                        "-h help");
                MPI_Finalize();
                return 0;
            default:
                fprintf(stderr, "Unknown option %s\n", optarg);
                break;
        }
    }


    if ((use_any_root && (root < 0 || root > num_proc)) || (!use_any_root && root)) {
        if (!rank)
            fprintf(stdout, "Wrong root\n");
        MPI_Finalize();
        return -1;
    }

    if (debug > 0) {
        fprintf(stdout, "%d: %s\n", rank, p_name);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (rank == root) {
        fprintf(stdout, "Inputs: root=%d, m=%d, M=%d, f=%d, rec=%d, reps=%d, a=%d, d=%d\n", root, message_min, message_max, msg_factor, rec, reps, alg, debug);
        fprintf(stdout, "root  p   g    msg     time     rec   reps   rec_in  rec_out  rec_world \n");
    }

#if HAVE_MPIX_H
    if (debug == 2)
       MPIX_Get_property(comm_world, MPIDO_RECT_COMM, &rec_world);
#endif

    for (msg_size = message_min; msg_size < message_max + 1; msg_size *= msg_factor) {
        array = create_rand_elms(msg_size);

        elapsed_time = 0.;
        for (i = 0; i < reps; i++) {
            MPI_Barrier(MPI_COMM_WORLD);
            start_time = MPI_Wtime();
            hpnla_bcast(array, msg_size, MPI_CHAR, root, MPI_COMM_WORLD, alg);
            elapsed_time += (MPI_Wtime() - start_time);
        }
        elapsed_time /= reps; // time in sec

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);

        if (rank == root) {
            fprintf(stdout, "MPI_Bcast: %d %d %d %f %f %d %d %d %d %d\n",
                    root, num_proc, 1, msg_size * sizeof (char) / 1024., max_time, rec, reps, -1, -1, rec_world);
        }

        free(array);
    }


    if (rec != 0) {
        for (msg_size = message_min; msg_size < message_max + 1; msg_size *= msg_factor) {

			array = create_rand_elms(msg_size);
			MPI_Barrier(MPI_COMM_WORLD);

			elapsed_time = 0.;
			for (i = 0; i < reps; i++) {
				MPI_Barrier(MPI_COMM_WORLD);
				start_time = MPI_Wtime();
				MPI_HBcast(array, msg_size, MPI_CHAR, root, MPI_COMM_WORLD, rec, alg);
				elapsed_time += (MPI_Wtime() - start_time);
			}

			elapsed_time /= reps;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);

			if (rank == root) {
				fprintf(stdout, "MPI_HBcast: %d %d %d %f %f %d %d %d\n",
						root,
						num_proc,
						hmpi_get_comm_conf(MPI_COMM_WORLD, "fayil.conf")->num_groups,
						msg_size * sizeof (char) / 1024.,
						max_time,
						rec,
						reps,
						rec_world);
			}

			free(array);
			MPI_Barrier(MPI_COMM_WORLD);

        }
    }

    if (debug == 2) {        // benchmark mpi_comm_split
        MPI_Comm in_group_comm, out_group_comm;
        double max_time_out, max_time_in;

        if (num_proc > HBCAST_MIN_PROCS) {
            if (rank == root)
                fprintf(stdout, "SPLIT. root p g max_time_in max_time_out rec reps\n");

            groups = 2;
            while (groups < g_max) {
                if (validate_groups(groups, num_proc) && num_proc % groups == 0) {
                    pg = num_proc / groups;
                    my_group = rank / pg;

                    elapsed_time = 0.;
                    for (i = 0; i < reps; i++) {
                        MPI_Barrier(MPI_COMM_WORLD);
                        start_time = MPI_Wtime();
                        MPI_Comm_split(MPI_COMM_WORLD, (rank - my_group * pg == 0) ? 0 : MPI_UNDEFINED, rank, &out_group_comm);
                        elapsed_time += (MPI_Wtime() - start_time);
                    }
                    elapsed_time /= reps;
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Reduce(&elapsed_time, &max_time_out, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);

                    elapsed_time = 0.;
                    for (i = 0; i < reps; i++) {
                        MPI_Barrier(MPI_COMM_WORLD);
                        start_time = MPI_Wtime();
                        MPI_Comm_split(MPI_COMM_WORLD, my_group, rank, &in_group_comm);
                        elapsed_time += (MPI_Wtime() - start_time);
                    }
                    elapsed_time /= reps;
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Reduce(&elapsed_time, &max_time_in, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);

                    if (rank == root) {
                        fprintf(stdout, "%d %d %d %f %f %d %d\n",
                                root,
                                num_proc,
                                groups,
                                max_time_in,
                                max_time_out,
                                rec,
                                reps
                                );
                    }

                    if (out_group_comm != MPI_COMM_NULL)
                        MPI_Comm_free(&out_group_comm);
                    if (in_group_comm != MPI_COMM_NULL)
                        MPI_Comm_free(&in_group_comm);
                }
                groups++;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    /* Shut down MPI */
    MPI_Finalize();
} /* main */
