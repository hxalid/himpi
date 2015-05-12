/*
Actual broadcast a-la-Torrent style
Author: Kiril Dichev
*/

#include "mpi.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#define CHUNK_SIZE (256 * 16)
#define CHUNKS 1024

pthread_mutex_t mutex_got_bit = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_ready = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t signal_got_bit = PTHREAD_COND_INITIALIZER;
pthread_cond_t signal_ready = PTHREAD_COND_INITIALIZER;
struct mpi_data {
int rank;
int size;
MPI_Comm comm;
double time1;
pthread_mutex_t *mutex_indexes;
};
typedef struct mpi_data mpi_data_t;
static pthread_t threads[2];
int indexes[CHUNKS];	
int alien_indexes[CHUNKS];	
char local_buffer[CHUNKS * CHUNK_SIZE];	

int done(int indexes[CHUNKS], pthread_mutex_t *mutex_indexes) {
	int i;
	int done = 1;
	pthread_mutex_lock(mutex_indexes);
	for (i=0; i < CHUNKS; i++) {
		if (indexes[i] == 0) {
			done = 0;
			break;
		}
	}
	pthread_mutex_unlock(mutex_indexes);
	return done;
} 

int has_interesting_bits(int *alien_indexes, int *indexes, pthread_mutex_t *mutex_indexes) {
	int i;
	int bits = -1;
	pthread_mutex_lock(mutex_indexes);
	for (i=0; i < CHUNKS; i++) {
		if (alien_indexes[i] > indexes[i]) {
			bits = i;
			break;
		}
	}
	pthread_mutex_unlock(mutex_indexes);
	return bits;
}


/*RESPONDER LOOP*/
void *signal_arrival_of_data_loop(void *args) {
	mpi_data_t *t = args;
	MPI_Status status;
	int i,j;
	int *copy_indexes;
	/*
	peers that want data from me
	*/
	int interested_peers[4];
	//the peers
	for (i=0;i<4;i++) {
		interested_peers[i] = (t->rank + t->size - i - 1) % t->size;
	}
	MPI_Request request[4];
	//if you are root, signal the interested peers that you have all bits of data
	if (t->rank == 0) {
		for (j=0; j<CHUNKS;j++) {
			for (i=0; i<4; i++) {
				MPI_Isend(indexes, CHUNKS, MPI_INT, interested_peers[i], 2, t->comm, &request[i]);
			}
			MPI_Waitall(4, request, MPI_STATUSES_IGNORE);
		}
	}
	//else, wait until data arrives in local buffer; when it does - signal the peers interested in you
	else {
		while (1) {

			//this mutex is held for long because otherwise there's multiple signals for 1 wait
			//which breaks the logic of the program
			pthread_mutex_lock(&mutex_got_bit);
			pthread_cond_wait(&signal_got_bit, &mutex_got_bit);
			copy_indexes = (int *) malloc(CHUNKS * sizeof(int));
			pthread_mutex_lock(t->mutex_indexes);
			for (j=0; j< CHUNKS; j++) 
				copy_indexes[j] = indexes[j];
			pthread_mutex_unlock(t->mutex_indexes);
			for (i=0;i<4;i++) {
				if (interested_peers[i] != 0) {
					MPI_Isend(copy_indexes, CHUNKS, MPI_INT, interested_peers[i], 2, t->comm, &request[i]);
				}
				else {
					request[i] = MPI_REQUEST_NULL;
				}
			}
			MPI_Waitall(4, request, MPI_STATUSES_IGNORE);
			free(copy_indexes); copy_indexes = NULL;
			pthread_mutex_unlock(&mutex_got_bit);
		}
	}
}

void *terminator_loop(void * args) {

	mpi_data_t *t = args;
	int flag = 0;
	if (t->rank != 0) {
		pthread_mutex_lock(&mutex_ready);
		pthread_cond_wait(&signal_ready, &mutex_ready);
		pthread_mutex_unlock(&mutex_ready);
	}
	MPI_Barrier(t->comm);
	double time2 = MPI_Wtime() - t->time1;
	double max_time;
	MPI_Reduce(&time2, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, t->comm);
	if (t->rank == 0) {
		printf("TOTAL TIME: %lf\n", max_time);
	}
	//go_on = 0;
	//signal the other thread on same node to stop
	MPI_Send(NULL, 0, MPI_INT, t->rank, 999, t->comm);
	
	pthread_exit(0);

}

/*REQUESTOR LOOP*/
void *requester_loop(void *args) {
	
	mpi_data_t *t = (mpi_data_t *) args;
	MPI_Status status;
	int i;
	int interesting_bit;
	int request=0;

	pthread_t thread1, thread2;
	pthread_create(&thread1, NULL, signal_arrival_of_data_loop, args);
	pthread_create(&thread2, NULL, terminator_loop, args);


	if (t->rank == 0) {
		while (1) {
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, t->comm, &status);
			if (status.MPI_TAG == 3) {
				MPI_Recv(&interesting_bit, 1, MPI_INT, MPI_ANY_SOURCE, 3, t->comm, &status);
				MPI_Send(&local_buffer[interesting_bit * CHUNK_SIZE], CHUNK_SIZE, MPI_CHAR, status.MPI_SOURCE, 4, t->comm);
			}
			else if ( status.MPI_TAG == 999) {
				MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, 999, t->comm, &status);
				pthread_exit(0);
			}
			else {
				printf("master got message with tag %d\n", status.MPI_TAG);
			}
		}
	}
	else {
		while  (1) {
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, t->comm, &status);
			if (status.MPI_TAG == 3) {
				MPI_Recv(&interesting_bit, 1, MPI_INT, MPI_ANY_SOURCE, 3, t->comm, &status);
				MPI_Send(&local_buffer[interesting_bit * CHUNK_SIZE], CHUNK_SIZE, MPI_CHAR, status.MPI_SOURCE, 4, t->comm);
			}
			else if ( status.MPI_TAG == 999) {
				MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, 999, t->comm, &status);
				pthread_exit(0);
			}
			else if (status.MPI_TAG == 2) {
				MPI_Recv(alien_indexes, CHUNKS, MPI_INT, MPI_ANY_SOURCE, 2, t->comm, &status);
				interesting_bit = has_interesting_bits(alien_indexes, indexes, t->mutex_indexes);
				if (interesting_bit != -1) {
					//request data from peer
					MPI_Sendrecv(&interesting_bit, 1, MPI_INT, status.MPI_SOURCE, 3, &local_buffer[interesting_bit *  CHUNK_SIZE],  CHUNK_SIZE , MPI_CHAR, status.MPI_SOURCE, 4, t->comm, MPI_STATUS_IGNORE);
					//record that bit was received
					pthread_mutex_lock(t->mutex_indexes);
					indexes[interesting_bit] = 1;
					pthread_mutex_unlock(t->mutex_indexes);
					//print(local_buffer);
					//BUG: IF CALLED TWICE FOR EACH CORRESPONDING COND_WAIT, ONE WAIT WILL BE PASSED
					pthread_mutex_lock(&mutex_got_bit);
					if (0 != pthread_cond_signal( &signal_got_bit)) printf("SIGNAL DID NOT RETURN 0\n");
					pthread_mutex_unlock(&mutex_got_bit);

					//signal to thread which handles 4 interested peers (not 4 sender peers)
				}
			}
			else {
				printf("got message with TAG = %d\n", status.MPI_TAG);
			}
			if (done(indexes, t->mutex_indexes)) {
				pthread_mutex_lock(&mutex_ready);
				if (0 != pthread_cond_signal( &signal_ready)) printf("SIGNAL DID NOT RETURN 0\n");
				pthread_mutex_unlock(&mutex_ready);
			}
		}
	}
}


void print(char *local_buffer) {
	int i;
	for (i=0; i < (CHUNKS * CHUNK_SIZE); i++) {
		printf( "%c", local_buffer[i]);
	}
	printf("\n");
	fflush(stdout);
}

int main(int argc, char *argv[]) {
	//MPI_Request reqs[4];
	int rank, size, i, provided;
	void *status;
	int required = MPI_THREAD_MULTIPLE;
	pthread_mutex_t *mutex_indexes;
	MPI_Init_thread(&argc, &argv, required, &provided);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	mpi_data_t mpi_data;
	mpi_data.rank = rank;
	mpi_data.size = size;
	mpi_data.comm = MPI_COMM_WORLD;
	//mutex_indexes
	mutex_indexes = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t));
	pthread_mutex_init(mutex_indexes, NULL);
	mpi_data.mutex_indexes = mutex_indexes;

	if (provided !=  required) {
		printf("Provided = %d,required = %d \n", provided, required);
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	//Initialize data and indexes to be broadcast
	if (rank == 0) {
		for (i=0;i < CHUNKS; i++) indexes[i] = 1;
		for (i=0;i< CHUNKS * CHUNK_SIZE;i++) local_buffer[i] = '0';
	}
	else {
		for (i=0;i < CHUNKS; i++) indexes[i] = 0;
		for (i=0;i< CHUNKS * CHUNK_SIZE;i++) local_buffer[i] = '1';
	}

	mpi_data.time1 = MPI_Wtime();

	pthread_create(&threads[0], NULL, requester_loop, &mpi_data);
	pthread_join(threads[0], &status); 

	pthread_mutex_destroy(mutex_indexes);

	printf("before finalize ...\n");
	MPI_Finalize();
}
