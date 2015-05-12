/*
 * p2p-pairs.c
 *
 *  Created on: 25 Aug 2011
 *      Author: kiril
 */

#include "config.h"

#include <getopt.h>
//#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "benchmarks/mpib.h"
#include "p2p/mpib_p2p.h"

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	//MPIB_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;

    //try to read a hosts file which explicitly gives the mapping process-host
    FILE *fp;
    if (NULL != (fp = fopen("hosts", "r"))) {
            MPI_Comm newcomm;
            char name[MPI_MAX_PROCESSOR_NAME];
            int len;
            MPI_Get_processor_name(name, &len);
            int color = 1;
            int key = 0;

            char line[MPI_MAX_PROCESSOR_NAME];
            while (fgets(line,sizeof(line),fp)) {
                    if (strncmp(line, name,len) == 0) {
                            break;
                    }
                    else
                            key++;
            }
            MPI_Comm_split(comm, color, key, &newcomm);
            comm = newcomm;
    }

	int exit;
	MPIB_getopt_help_default(&exit);
	MPIB_msgset msgset;
	MPIB_getopt_msgset_default(&msgset);
	// if we want adaptive message sizes, we should do: msgset.stride = 0;
	MPIB_precision precision;
	MPIB_getopt_precision_default(&precision);
	int rank;
	MPI_Comm_rank(comm, &rank);

	if (rank == 0) {
		int c;
		char options[256] = "i:P:";
		strcat(options, MPIB_getopt_help_options());
		strcat(options, MPIB_getopt_msgset_options());
		strcat(options, MPIB_getopt_precision_options());
		strcat(options, MPIB_getopt_p2p_options());
		while ((c = getopt(argc, argv, options)) != -1) {
			MPIB_getopt_help_optarg(c, "p2p", &exit);
			MPIB_getopt_msgset_optarg(c, &msgset);
			MPIB_getopt_precision_optarg(c, &precision);
		}
	}
	MPIB_getopt_help_bcast(&exit, 0, comm);
	if (exit) {
		MPI_Finalize();
		return 0;
	}
	MPIB_getopt_msgset_bcast(&msgset, 0, comm);
	MPIB_getopt_precision_bcast(&precision, 0, comm);


	int size;
	MPI_Comm_size(comm, &size);
	if ((size % 2) != 0) {
		if (rank == 0)
			fprintf(stderr, "This tool is currently tested only for an even number of procs!\n");
		MPI_Abort(comm, -1);
	}
	else {
		if (rank == 0)
			fprintf(stdout,
"# The pairs should be started with the even numbers at one partition, the odd numbers at the other partition\n");
	}
	//Since in the two pairs there are two senders, we halve the message size
	//This way, the total transferred message size is as specified by the user
	msgset.min_size = (2. * msgset.min_size ) / (double) size;
	msgset.max_size = (2. * msgset.max_size) / (double) size;
	msgset.stride = (2. * msgset.stride) / (double) size;


	if (rank == 0)
		printf("#p2p adaptive benchmark 0-1\n#\n");
	MPIB_print_processors(comm);
	if (rank == 0) {
		MPIB_print_msgset(msgset);
		MPIB_print_precision(precision);
	}

	int count = 0;
	MPIB_result* results = NULL;
	double time = MPI_Wtime();

	MPIB_measure_p2p_parallel_msgset(MPI_Send, MPI_Recv, comm, msgset, precision, &count, &results);

	if (rank == 0) {
		printf("#total benchmarking time: %le\n#\n", MPI_Wtime() - time);
		MPIB_print_result_th();
		int i;
		for (i = 0; i < count; i++) {
			MPIB_print_result_tr(results[i]);
		}
	}
	free(results);

	MPI_Finalize();
	return 0;
}

