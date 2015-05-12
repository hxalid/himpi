#include "mpib_p2p_benchmarks.h"
#include <stdio.h>
//#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "mpib_output.h"
#include "mpib_utilities.h"

int MPIB_check_comm_size(MPI_Comm comm) {
	int rank;
	MPI_Comm_rank(comm, &rank);
	int size;
	MPI_Comm_size(comm, &size);
	if (size < 2) {
		if (rank == 0)
			fprintf(stderr, "Cannot perform p2p benchmarks for %d processes (must be >= 2)\n", size);
		return -1;
	}
	return 0;
}

void MPIB_measure_p2p(MPIB_Send send, MPIB_Recv recv, MPI_Comm comm, int measure, int mirror,
		int M, MPIB_precision precision, MPIB_result* result) {
	if (MPIB_check_comm_size(comm))
		return;
	int rank;
	MPI_Comm_rank(comm, &rank);
	if (rank == measure || rank == mirror) {
		double* T = rank == measure ? (double*)malloc(sizeof(double) * precision.max_reps) : NULL;
		char* buffer = (char*)malloc(sizeof(char) * M);
		int stop = 0;
		double sum = 0;
		int reps = 0;
		double ci = 0;
		while (!stop && reps < precision.max_reps) {
			if (rank == measure) {
				MPI_Recv(NULL, 0, MPI_CHAR, mirror, 0, comm, MPI_STATUS_IGNORE);
				T[reps] = MPI_Wtime();
				send(buffer, M, MPI_CHAR, mirror, 0, comm);
				recv(buffer, M, MPI_CHAR, mirror, 0, comm, MPI_STATUS_IGNORE);
				sum += T[reps] = MPI_Wtime() - T[reps];
			}
			else {
				MPI_Send(NULL, 0, MPI_CHAR, measure, 0, comm);
				recv(buffer, M, MPI_CHAR, measure, 0, comm, MPI_STATUS_IGNORE);
				send(buffer, M, MPI_CHAR, measure, 0, comm);
			}
			reps++;
			if (reps >= precision.min_reps && reps > 2) {
				if (rank == measure) {
					stop = (ci = MPIB_ci(precision.cl, reps, T)) * reps / sum < precision.eps;
					MPI_Send(&stop, 1, MPI_INT, mirror, 0, comm);
				}
				else
					MPI_Recv(&stop, 1, MPI_INT, measure, 0, comm, MPI_STATUS_IGNORE);
			}
		}
		free(buffer);
		if (rank == measure)
			*result = (MPIB_result){M, sum / reps, MPI_Wtick(), reps, ci};
	}
}

void MPIB_measure_p2p_msgset(MPIB_Send send, MPIB_Recv recv, MPI_Comm comm, int measure, int mirror,
		MPIB_msgset msgset, MPIB_precision precision, int* count, MPIB_result** results)
{
	if (MPIB_check_comm_size(comm))
		return;
	int rank;
	MPI_Comm_rank(comm, &rank);
	if (rank == measure || rank == mirror) {
		if (rank == measure) {
			if (MPIB_verbose)
				fprintf(stderr, "adaptive p2p benchmark %d-%d:", measure, mirror);
			*count = 0;
			*results = NULL;
		}
		int M = msgset.min_size;
		int c = 0;
		// variables at measure
		double* T = NULL;
		double wtick = 0;
		int stride = msgset.min_stride; //adaptive stride
		int pos = 0;
		if (rank == measure) {
			// use the same array and wtick value for communication experiments with different message sizes
			T = (double*)malloc(sizeof(double) * precision.max_reps);
			wtick = MPI_Wtick();
		}
		/*
		In any case, we don't do more than max_size message sizes, but
		with the adaptive approach we also don't iterate more than max_num times
		*/
		while (M <= msgset.max_size && !( (msgset.stride == 0) && (c < msgset.max_num))) {
			if ((rank == measure) && MPIB_verbose) {
				fprintf(stderr, " %d", M);
				fflush(stderr);
			}
			char* buffer = (char*)malloc(sizeof(char) * M);
			int stop = 0;
			double sum = 0;
			int reps = 0;
			double ci = 0;
			/*
			This loop performs a number of repetitions with a fixed M. It runs while
			max(min_reps,2) < reps < max_reps) and NOT (ci * reps / sum ) < eps
			*/
			while (!stop && reps < precision.max_reps) {
				if (rank == measure) {
					MPI_Recv(NULL, 0, MPI_CHAR, mirror, 0, comm, MPI_STATUS_IGNORE);
					T[reps] = MPI_Wtime();
					send(buffer, M, MPI_CHAR, mirror, 0, comm);
					recv(buffer, M, MPI_CHAR, mirror, 0, comm, MPI_STATUS_IGNORE);
					sum += T[reps] = MPI_Wtime() - T[reps];
				}
				else {
					MPI_Send(NULL, 0, MPI_CHAR, measure, 0, comm);
					recv(buffer, M, MPI_CHAR, measure, 0, comm, MPI_STATUS_IGNORE);
					send(buffer, M, MPI_CHAR, measure, 0, comm);
				}
				reps++;
				if (reps >= precision.min_reps && reps > 2) {
					if (rank == measure) {
						stop = (ci = MPIB_ci(precision.cl, reps, T)) * reps / sum < precision.eps;
						MPI_Send(&stop, 1, MPI_INT, mirror, 0, comm);
					}
					else
						MPI_Recv(&stop, 1, MPI_INT, measure, 0, comm, MPI_STATUS_IGNORE);
				}
			}
			free(buffer);
			c++;
			if (rank == measure) {
				MPIB_result result = (MPIB_result){M, sum / reps, wtick, reps, ci};
				*count = c;
				if (pos < c - 1) {
					MPIB_result* tmp = (MPIB_result*)malloc(sizeof(MPIB_result) * c);
					memcpy(tmp, *results, sizeof(MPIB_result) * pos);
					memcpy(&tmp[pos + 1], &(*results)[pos], sizeof(MPIB_result) * (c - 1 - pos));
					free(*results);
					*results = tmp;
				}
				else
					*results = (MPIB_result*)realloc(*results, sizeof(MPIB_result) * c);
				(*results)[pos] = result;
				/*
				determine the next M - adaptively or with fixed stride:
				if msgset.stride > 0, then use this fixed stride (non-adaptive)
				if msgset.stride == 0, increment M adaptively 
				*/
				if (msgset.stride == 0) {
					//adaptive M increment :
					do {
						if (pos > 1 && (MPIB_diff(result, &(*results)[pos - 2]) < msgset.max_diff)) {
							stride *= 2;
							M += stride;
							while (pos < c && (*results)[pos].M < M)
								pos++;
						} else {
							if ((stride / 2) > msgset.min_stride) {
								stride /= 2;
								M -= stride;
							} else {
								M += stride;
								while (pos < c && (*results)[pos].M < M)
									pos++;
							}
						}
					} while (pos < c && (*results)[pos].M == M);
				} else {
					M += msgset.stride;
					pos++;
				}
				MPI_Send(&M, 1, MPI_INT, mirror, 0, comm);
			} else
				MPI_Recv(&M, 1, MPI_INT, measure, 0, comm, MPI_STATUS_IGNORE);
		}
		free(T);
		if ((rank == measure) && MPIB_verbose)
			fprintf(stderr, "\n\n");
	}
}

void MPIB_measure_p2p_parallel_msgset(MPIB_Send send, MPIB_Recv recv, MPI_Comm comm,
		MPIB_msgset msgset, MPIB_precision precision, int* count, MPIB_result** results) {

	MPI_Comm comm_1;
	int key = MPI_UNDEFINED;
	int rank;
	MPI_Comm_rank(comm, &rank);
	if ((rank % 2) == 0) {
		key = 0;
	}
	//new_comm will be NULL for processes on receiver side and contain all processes on sender side for the senders
	MPI_Comm_split(comm, key, rank, &comm_1);

	// variables at measure
	double* T = NULL;
	double wtick = 0;
	if (comm_1 != MPI_COMM_NULL) {
		// use the same array and wtick value for communication experiments with different message sizes
		T = (double*)malloc(sizeof(double) * precision.max_reps);
		wtick = MPI_Wtick();
		*count = 0;
		*results = NULL;
	}

	int M = msgset.min_size;
	int c = 0;
	/*
	Don't do more than max_size message sizes. The adaptive approach is ignored in this routine due to the number
	of communicating pairs which need to run the same number of iterations
	 */
	while ((M <= msgset.max_size) && (c < msgset.max_num)) {
		char* buffer = (char*)malloc(sizeof(char) * M);
		int reps = 0;
		double time_over_reps;
		/*
		This loop performs a number of repetitions with a fixed M. It runs while
		max(min_reps,2) < reps < max_reps) and NOT (ci * reps / sum ) < eps
		 */
		for (reps = 0; reps < precision.max_reps; reps++) {
			if (comm_1 != MPI_COMM_NULL) {
				//all sending processes synchronize before sending data
				MPI_Barrier(comm_1);
				T[reps] = MPI_Wtime();
				send(buffer, M, MPI_CHAR, rank+1, 0, comm);
				recv(buffer, M, MPI_CHAR, rank+1, 0, comm, MPI_STATUS_IGNORE);
				T[reps] = MPI_Wtime() - T[reps];
				//send the maximum time of all sender processes to process 0
				if (rank == 0)
					MPI_Reduce(MPI_IN_PLACE, &T[reps], 1, MPI_DOUBLE, MPI_MAX, 0, comm_1);
				else
					MPI_Reduce(&T[reps], NULL, 1, MPI_DOUBLE, MPI_MAX, 0, comm_1);
			}
			else {
				recv(buffer, M, MPI_CHAR, rank-1, 0, comm, MPI_STATUS_IGNORE);
				send(buffer, M, MPI_CHAR, rank-1, 0, comm);
			}
		}
		free(buffer);
		c++;
		*count = c;
		if (rank == 0) {
			//strategy over repetitions - max, min or average (currently min)
			time_over_reps = DBL_MAX;
			for (reps = 0; reps < precision.max_reps; reps++) {
				time_over_reps = (T[reps] < time_over_reps)?T[reps]:time_over_reps;
			}
			*results = (MPIB_result*)realloc(*results, sizeof(MPIB_result) * c);
			int size;
			MPI_Comm_size(comm, &size);
			size /=2;
			MPIB_result result = (MPIB_result){M * size /*total message size between partitions*/, time_over_reps, wtick, reps, /*ci*/0.};
			(*results)[c-1] = result;
		}
		M += msgset.stride;
	}
	free(T);

}

void MPIB_measure_allp2p(MPIB_Send send, MPIB_Recv recv, MPI_Comm comm, int parallel,
		int M, MPIB_precision precision, MPIB_result* results)
{
	if (MPIB_check_comm_size(comm))
		return;
	int rank;
	MPI_Comm_rank(comm, &rank);
	int size;
	MPI_Comm_size(comm, &size);
	int sendcount =  size - rank - 1;
	// results for pairs including a given processor
	MPIB_result* res = (MPIB_result*)malloc(sizeof(MPIB_result) * sendcount);
	// use the same array and wtick value for communication experiments with different pairs including a given processor
	double* T = (double*)malloc(sizeof(double) * precision.max_reps);
	double wtick = MPI_Wtick();
	char* buffer = (char*)malloc(sizeof(char) * M);
	MPIB_pairs* pairs = MPIB_build_pairs(size);
	MPIB_pairs* iter_pairs = pairs;
	while (iter_pairs) {
		MPI_Barrier(comm);
		MPIB_pair* iter_pair = iter_pairs->list;
		while (iter_pair) {
			if (!parallel)
				MPI_Barrier(comm);
			int i = iter_pair->values[0];
			int j = iter_pair->values[1];
			if (rank == i || rank == j) {
				int stop = 0;
				double sum = 0;
				int reps = 0;
				double ci = 0;
				while (!stop && reps < precision.max_reps) {
					if (rank == i) {
						MPI_Recv(NULL, 0, MPI_CHAR, j, 0, comm, MPI_STATUS_IGNORE);
						T[reps] = MPI_Wtime();
						send(buffer, M, MPI_CHAR, j, 0, comm);
						recv(buffer, M, MPI_CHAR, j, 0, comm, MPI_STATUS_IGNORE);
						sum += T[reps] = MPI_Wtime() - T[reps];
					}
					else {
						MPI_Send(NULL, 0, MPI_CHAR, i, 0, comm);
						recv(buffer, M, MPI_CHAR, i, 0, comm, MPI_STATUS_IGNORE);
						send(buffer, M, MPI_CHAR, i, 0, comm);
					}
					reps++;
					if (reps >= precision.min_reps && reps > 2) {
						if (rank == i) {
							stop = (ci = MPIB_ci(precision.cl, reps, T)) * reps / sum < precision.eps;
							MPI_Send(&stop, 1, MPI_INT, j, 0, comm);
						}
						else
							MPI_Recv(&stop, 1, MPI_INT, i, 0, comm, MPI_STATUS_IGNORE);
					}
				}
				if (rank == i)
					res[j - i - 1] = (MPIB_result){M, sum / reps, wtick, reps, ci};
			}
			iter_pair = iter_pair->next;
		}
		iter_pairs = iter_pairs->next;
	}
	MPIB_free_pairs(pairs);
	free(buffer);
	free(T);
	// allgather the results
	int* recvcounts = (int*)malloc(sizeof(int) * size);
	int* displs = (int*)malloc(sizeof(int) * size);
	int i;
	for (i = 0; i < size; i++) {
		recvcounts[i] = sizeof(MPIB_result) * (size - i - 1);
		displs[i] = i == 0 ? 0 : displs[i - 1] + recvcounts[i - 1]; 
	}
	MPI_Allgatherv(res, sizeof(MPIB_result) * sendcount, MPI_CHAR, results, recvcounts, displs, MPI_CHAR, comm);
	free(recvcounts);
	free(displs);
	free(res);
}
