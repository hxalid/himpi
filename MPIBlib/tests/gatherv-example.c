#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <limits.h>
#include <float.h>

#include "mpi.h"


#define BUF_SIZE 16384
#define NUM_FRAGMENTS 600
#define NUM_ITER 5


int main(int argc,char *argv[])
{
	int    n, myid, num_procs, i, j;
	double startwtime = 0.0, endwtime;
	int    namelen;
	FILE *fp = NULL;
	char *root_buffer;
	double root_totals[NUM_ITER];
	char   processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Comm comm;
	MPI_Init(&argc,&argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &num_procs);
	int my_fragments[num_procs];
	MPI_Comm_rank(comm, &myid);

	int total_fragments_num = 0;
	double factors[num_procs];

	//my_fragments initialization
	for (i=0; i<num_procs; i++) {
		my_fragments[i] = NUM_FRAGMENTS ;
	}
	//exclude own fragments, they are not used
	for (i=1; i< num_procs; i++) total_fragments_num += my_fragments[i]; 

	if (myid == 0) {
		for (i = 0; i < num_procs; i++) factors[i] = 1;

		int c;
		char options[256] = "f:";
		while ((c = getopt(argc, argv, options)) != -1) {
			switch (c) {
				case 'f': 
					fp = fopen(optarg, "r");
					if (fp == NULL) {
						fprintf(stderr, "Cannot open file %s for reading\n", optarg);
						MPI_Abort(comm, -1);
					} else {
						int i;
						for (i = 0; i < num_procs; i++)
							fscanf(fp, "%le\n", &factors[i]);
						fclose(fp);
					}
			}
		}
	}
	MPI_Bcast(factors, num_procs, MPI_DOUBLE, 0, comm);

	if (myid == 0) {
		root_buffer = (char *) malloc (sizeof(char) * total_fragments_num * BUF_SIZE);
	}

	char buffer[BUF_SIZE];
	for (j=0; j < NUM_ITER; j++) {
		for (i=0; i< BUF_SIZE; i++) buffer[i] = 's';
		
		MPI_Barrier(comm);
		double start = MPI_Wtime();
		//one gatherv operation
		if (myid != 0) {
			MPI_Request req[my_fragments[myid]];
			for (i=0;i<my_fragments[myid];i++) {
				MPI_Isend(buffer, BUF_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &req[i]);
			}
			MPI_Waitall(my_fragments[myid], req, MPI_STATUSES_IGNORE);
		}
		else {
			MPI_Request req[total_fragments_num];
			char *ptr = root_buffer;
			for (i=0;i<total_fragments_num;i++) {
				MPI_Irecv(ptr, BUF_SIZE, MPI_CHAR, MPI_ANY_SOURCE, 0, comm, &req[i]);
				ptr += BUF_SIZE;
			}
			//get all fragments from all peers in any order
			MPI_Waitall(total_fragments_num, req, MPI_STATUSES_IGNORE);
		}

		double end = MPI_Wtime();
		double total = end - start;
		double root_total;
		MPI_Reduce(&total, &root_total, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
			root_totals[j] = root_total;
	}
	if (myid == 0) {
		double min_of_maxes = DBL_MAX;
		for (j=0;j<NUM_ITER;j++) {
			//printf("root_totals[%d] = %lf\n", j, root_totals[j]);
			if (root_totals[j] < min_of_maxes)
				min_of_maxes = root_totals[j];
		}
		printf("time = %lf\n", min_of_maxes);
	}
	if (myid == 0) {
		free(root_buffer);
		root_buffer = NULL;
	}
	MPI_Finalize();
	return 0;
}
