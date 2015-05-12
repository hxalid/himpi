/*!
 * \page p2p P2P benchmark
 * Performs adaptive p2p benchmark between a pair of processors 0-1.
 * Arguments are described in \ref getopt.
 * 
 * \b Example 1 (standard P2P benchmark):
 * \code
 * $ mpirun -n 2 p2p -m 1024 -M 2049 > p2p.out
 * \endcode
 * will benchmark the p2p communication between 2 processes for
 * message sizes of 1024 and 2048 bytes
 *
 * \b Example 2 (scatter-gather-based p2p benchmarks, mainly useful for
 * misconfigured long-haul connections with high-bandwidth and high-latency):
 * \code
 * cd <install-dir>/bin
 * mpirun -n 2 p2p -m 1024 -M 2029 -P 1 > p2p.out
 * \endcode
 * This example will run with default of 2 extra processes at sender
 * and receiver side, respectively.
 * 
 * \b p2p.plot draws the graph of the execution time (sec) against message size (kb).
 * - input: p2p.out
 * - output: p2p.eps (0-1 with error bars)
 *
 * Using the gnuplot script:
 * \code
 * $ gnuplot p2p.plot
 * \endcode
 */

#include "config.h"

#include <getopt.h>
//#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "benchmarks/mpib.h"
#include "p2p/mpib_p2p.h"

extern MPI_Comm intracomm;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	//MPIB_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;

	int exit;
	int p2p_type;
	int parallel = 0;
	int proc_num_spawned = 0;
	MPIB_getopt_help_default(&exit);
	MPIB_getopt_p2p_default(&parallel, &p2p_type, &proc_num_spawned);
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
			MPIB_getopt_p2p_optarg(c, &parallel, &p2p_type, &proc_num_spawned);
		}
	}
	MPIB_getopt_help_bcast(&exit, 0, comm);
	if (exit) {
		MPI_Finalize();
		return 0;
	}
	MPIB_getopt_msgset_bcast(&msgset, 0, comm);
	MPIB_getopt_precision_bcast(&precision, 0, comm);
	MPIB_getopt_p2p_bcast(&parallel, &p2p_type, &proc_num_spawned, 0, comm);

	if (parallel) {
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
	}

	if (rank == 0)
		printf("#p2p adaptive benchmark 0-1\n#\n");
	MPIB_print_processors(comm);
	if (rank == 0) {
		MPIB_print_msgset(msgset);
		MPIB_print_precision(precision);
	}

	p2p_init_max_fragment_size();

	if (p2p_type == 1) {
		p2p_init(comm, proc_num_spawned);
		char name[128];
		int real_length;
		int rank;
		MPI_Get_processor_name(name, &real_length);
		if (intracomm != MPI_COMM_NULL) {
			MPI_Comm_rank(intracomm, &rank);
			printf("I am an original process, rank is %d, name is %s\n", rank, name);
			fflush(stdout);
		}
		else
			printf("I am an original process, I don't have a rank, name is %s\n", name);
			
	}
	int count = 0;
	MPIB_result* results = NULL;
	double time = MPI_Wtime();

	MPIB_measure_p2p_msgset(MPI_Send, MPI_Recv, comm, 0, 1, msgset, precision, &count, &results);



	if (rank == 0) {
		printf("#total benchmarking time: %le\n#\n", MPI_Wtime() - time);
		MPIB_print_result_th();
		int i;
		for (i = 0; i < count; i++) {
			MPIB_print_result_tr(results[i]);
		}
	}
	free(results);
	
	if (p2p_type == 1) 
		p2p_finalize();
	MPI_Finalize();
	return 0;
}
