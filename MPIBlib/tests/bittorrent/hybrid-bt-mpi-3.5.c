/*
Actual broadcast a-la-Torrent style
Author: Kiril Dichev
*/

#include "mpi.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#define CHUNK_SIZE 32768
#define CHUNKS 5096
#define NUM_CHANNELS 1

pthread_mutex_t mutex_got_bit = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_indexes = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_ready = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t signal_ready = PTHREAD_COND_INITIALIZER;
struct mpi_data {
int rank;
int size;
MPI_Comm comm;
double time1;
};
typedef struct mpi_data mpi_data_t;
static pthread_t threads[2];
int indexes[CHUNKS];	
int alien_indexes[NUM_CHANNELS][CHUNKS];	
char local_buffer[CHUNKS * CHUNK_SIZE];	
int interested_peers[NUM_CHANNELS];
int interesting_peers[NUM_CHANNELS]; //as opposed to interested_peers
int last_received_index;
static double recv_timer = 0.;
static double send_timer = 0.;

void shuffle(int *array, size_t n) {
	srand(time(NULL));
	if (n > 1) {
		size_t i;
		for (i = 0; i < n - 1; i++) {
			size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
			int t = array[j];
			array[j] = array[i];
			array[i] = t;
		}
	}
}

void master_send_shuffled(struct mpi_data *t) {
	int shuffled_indexes[CHUNKS*NUM_CHANNELS];
	int i,j;
	for (j=0; j < CHUNKS*NUM_CHANNELS; j++) {
		shuffled_indexes[j] = (j % CHUNKS);
	}
	shuffle(shuffled_indexes, CHUNKS*NUM_CHANNELS);

	MPI_Request request[CHUNKS*NUM_CHANNELS];
	int k = 0;
	for (i=0; i<CHUNKS;i++) {
		for (j=0; j<NUM_CHANNELS; j++) {
			if (interested_peers[j] != t->rank) {
				MPI_Isend(shuffled_indexes+k, 1, MPI_INT, interested_peers[j], 2, t->comm, &request[k]);
				k++;
				//MPI_Send(shuffled_indexes+i, 1, MPI_INT, interested_peers[j], 2, t->comm);
			}
			else {
				request[k] = MPI_REQUEST_NULL;
				k++;
			}
		}
	}
	MPI_Waitall(CHUNKS*NUM_CHANNELS, request, MPI_STATUSES_IGNORE);
}

int done() {
	int i;
	int flag = 1;
	pthread_mutex_lock(&mutex_indexes);
	for (i=0; i < CHUNKS; i++) {
		if (indexes[i] == 0) {
			flag = 0;
			break;
		}
	}
	pthread_mutex_unlock(&mutex_indexes);
	return flag;
} 

int map(int rank, int size) {
	int i;
	int factor = 1;
	int a;
	for (i=0; i<NUM_CHANNELS; i++) {
		a = (rank - factor) % size;
		//interesting_peers[i] =  (a < 0) ? (a + size) : a;
		//factor = 2 * factor;
		interesting_peers[i] = (rank + (i + 1)) % size;
	}
	return 0;
}

int map2(int rank, int size) {
	int i;
	//the peers
	for (i=0;i<NUM_CHANNELS;i++) {
		interested_peers[i] = (rank + size - i - 1) % size;
		//interested_peers[i] = (t->rank + t->size - i - 8) % t->size;
		//interested_peers[i] = (t->rank + factor) % t->size;
		//factor = 2 * factor;
	}
	return 0;
}

int get_index_from_rank(int rank) {
	int i;
	for (i=0; i<NUM_CHANNELS; i++) {
		if (rank == interesting_peers[i]) return i;
	}
	printf("error at rank %d\n", rank);
	exit(-1);
}

int has_interesting_bits(int alien_index, int peer_rank) {
	int i;
	int bits = 0;
	pthread_mutex_lock(&mutex_indexes);
	int index = get_index_from_rank(peer_rank);
	if (!alien_indexes[index][alien_index]) alien_indexes[index][alien_index] = 1;
	if (alien_indexes[index][alien_index] > indexes[alien_index]) 
			bits = 1;
	pthread_mutex_unlock(&mutex_indexes);
	return bits;
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
	MPI_Reduce(&recv_timer, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, t->comm);
	if (t->rank == 0) {
		printf("TOTAL RECEIVE TIME: %lf\n", max_time);
	}
	MPI_Reduce(&send_timer, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, t->comm);
	if (t->rank == 0) {
		printf("TOTAL SEND TIME: %lf\n", max_time);
	}
	//go_on = 0;
	//signal the other thread on same node to stop
	MPI_Send(NULL, 0, MPI_INT, t->rank, 999, t->comm);

}


void *get_data(void *args) {

	MPI_Status status;
	mpi_data_t *t = (mpi_data_t *) args;
	int alien_index;
	MPI_Request req;
	
	while (1) {
		MPI_Recv(&alien_index, 1, MPI_INT, MPI_ANY_SOURCE, 2, t->comm, &status);
		int interesting_bit = has_interesting_bits(alien_index, status.MPI_SOURCE);
		if (interesting_bit != 0) {
			//request data from peer
			MPI_Isend(&alien_index, 1, MPI_INT, status.MPI_SOURCE, 3, t->comm, &req); 
			//MPI_Send(&alien_index, 1, MPI_INT, status.MPI_SOURCE, 3, t->comm); 
			double begin = MPI_Wtime();
			MPI_Recv(&local_buffer[alien_index *  CHUNK_SIZE],  CHUNK_SIZE , MPI_CHAR, status.MPI_SOURCE, 4, t->comm, MPI_STATUS_IGNORE);
			double end = MPI_Wtime();
			recv_timer += (end-begin);
			//record that bit was received
			pthread_mutex_lock(&mutex_indexes);
			indexes[alien_index] = 1;
			last_received_index = alien_index;
			pthread_mutex_unlock(&mutex_indexes);

			MPI_Request request[NUM_CHANNELS];
			int i;
			for (i=0;i<NUM_CHANNELS;i++) {
				if ((interested_peers[i] != 0) && (interested_peers[i] != t->rank)) {
					MPI_Isend(&last_received_index, 1, MPI_INT, interested_peers[i], 2, t->comm, &request[i]);
				}
				else {
					request[i] = MPI_REQUEST_NULL;
				}
			}
			MPI_Waitall(NUM_CHANNELS, request, MPI_STATUSES_IGNORE);

			if (done()) {
				pthread_mutex_lock(&mutex_ready);
				if (0 != pthread_cond_signal( &signal_ready)) printf("SIGNAL DID NOT RETURN 0\n");
				pthread_mutex_unlock(&mutex_ready);
			}
		}
	}
}
void *send_data(void *args) {
	mpi_data_t *t = (mpi_data_t *) args;
	MPI_Status status;
	int alien_index;
	MPI_Request req;
	while (1) {
		MPI_Recv(&alien_index, 1, MPI_INT, MPI_ANY_SOURCE, 3, t->comm, &status);
		double begin = MPI_Wtime();
		MPI_Isend(&local_buffer[alien_index * CHUNK_SIZE], CHUNK_SIZE, MPI_CHAR, status.MPI_SOURCE, 4, t->comm, &req);
		//MPI_Send(&local_buffer[alien_index * CHUNK_SIZE], CHUNK_SIZE, MPI_CHAR, status.MPI_SOURCE, 4, t->comm);
		double end = MPI_Wtime();
		send_timer += (end - begin);
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
	MPI_Init_thread(&argc, &argv, required, &provided);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	mpi_data_t mpi_data;
	mpi_data.rank = rank;
	mpi_data.size = size;
	mpi_data.comm = MPI_COMM_WORLD;
	pthread_t send_thread, get_thread, terminator_thread;

	if (provided !=  required) {
		printf("Provided = %d,required = %d \n", provided, required);
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	//Initialize data and indexes to be broadcast
	int j;
	for (i=0; i<NUM_CHANNELS; i++){
		for (j=0; j<CHUNKS; j++)
			alien_indexes[i][j] = 0;
	}
	if (rank == 0) {
		for (i=0;i < CHUNKS; i++) indexes[i] = 1;
		for (i=0;i< CHUNKS * CHUNK_SIZE;i++) local_buffer[i] = '0';
	}
	else {
		for (i=0;i < CHUNKS; i++) indexes[i] = 0;
		for (i=0;i< CHUNKS * CHUNK_SIZE;i++) local_buffer[i] = '1';
	}

	mpi_data.time1 = MPI_Wtime();

	map(rank, size);
	map2(rank, size);

	if (rank == 0) {
		master_send_shuffled(&mpi_data);
	}
	pthread_create(&send_thread, NULL, send_data, &mpi_data);
	pthread_create(&get_thread, NULL, get_data, &mpi_data);
	pthread_create(&terminator_thread, NULL, terminator_loop, &mpi_data);
	MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, 999, mpi_data.comm, MPI_STATUS_IGNORE);
	//printf("process %d to cancel threads\n", mpi_data.rank);
	pthread_cancel(send_thread);
	pthread_cancel(get_thread);

	MPI_Finalize();
}