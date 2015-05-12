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
double time1;
typedef struct mpi_data mpi_data_t;
static pthread_t threads[2];
int indexes[CHUNKS];	
int alien_indexes[NUM_CHANNELS][CHUNKS];	
char local_buffer[CHUNKS * CHUNK_SIZE];	
int interested_peers[NUM_CHANNELS];
int interesting_peers[NUM_CHANNELS]; //as opposed to interested_peers
static double local_recv_time = 0.;
static double remote_recv_time = 0.;
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

void master_send_shuffled(MPI_Comm comm) {
	int shuffled_indexes[CHUNKS*NUM_CHANNELS];
	int rank;
	MPI_Comm_rank(comm, &rank);
	int i,j;
	for (j=0; j < CHUNKS*NUM_CHANNELS; j++) {
		shuffled_indexes[j] = (j % CHUNKS);
	}
	shuffle(shuffled_indexes, CHUNKS*NUM_CHANNELS);

	MPI_Request request[CHUNKS*NUM_CHANNELS];
	int k = 0;
	for (i=0; i<CHUNKS;i++) {
		for (j=0; j<NUM_CHANNELS; j++) {
			if ((interested_peers[j] != 0) && (interested_peers[j] != 1)) {
				MPI_Isend(shuffled_indexes+k, 1, MPI_INT, interested_peers[j], 2, comm, &request[k]);
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

int receiver (int rank) { 
if (rank % 2 == 0) 
	return (rank+1);
else 
	return rank;
}

//rank is odd, interesting peers are even and not (rank - 1)
int map(int rank, int size) {
	int i;
	//int factor = 1;
	//int a;
	for (i=0; i<NUM_CHANNELS; i++) {
		//a = (rank - factor) % size;
		//interesting_peers[i] =  (a < 0) ? (a + size) : a;
		//factor = 2 * factor;
		interesting_peers[i] = (2*size - 2*(i+1) + (rank - 1)) % size;
	}
	return 0;
}

//rank is even, interested peers are odd and not (rank + 1)
int map2(int rank, int size) {
	int i;
	//the peers
	for (i=0;i<NUM_CHANNELS;i++) {
		interested_peers[i] = (rank + 2 * (i+1) + 1) % size;
		//printf("rank %d; interested_peers[%d] = %d\n", rank, i, interested_peers[i]);
		//interested_peers[i] = (rank + size - i - 1) % size;
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
	if (alien_indexes[index][alien_index] > indexes[alien_index]) {
			indexes[alien_index] = 1; //not true yet, but I will get data soon and so will the even node on my node
			bits = 1;
	}
	pthread_mutex_unlock(&mutex_indexes);
	return bits;
}

/* even rank gets data from odd rank on the same node */
void *even_get_local_data(void *args) {

	MPI_Comm *comm = (MPI_Comm *) args;

	int rank;
	int alien_index;
	MPI_Comm_rank(*comm, &rank);
	while (1) {
		if (done()) {
			//printf("DEBUG: rank %d signals done\n", rank);
			pthread_mutex_lock(&mutex_ready);
			if (0 != pthread_cond_signal( &signal_ready)) printf("SIGNAL DID NOT RETURN 0\n");
			pthread_mutex_unlock(&mutex_ready);
		}

		double begin = MPI_Wtime();
		MPI_Recv(&alien_index, 1, MPI_INT, (rank+1), 5, *comm, MPI_STATUS_IGNORE);
		MPI_Recv(&local_buffer[alien_index *  CHUNK_SIZE],  CHUNK_SIZE , MPI_CHAR, (rank+1), 4, *comm, MPI_STATUS_IGNORE);
		//record that bit was received
		pthread_mutex_lock(&mutex_indexes);
		indexes[alien_index] = 1;
		pthread_mutex_unlock(&mutex_indexes);
		
		printf("DEBUG: got local data : %d , bit %d\n", rank, alien_index);
		double end = MPI_Wtime();
		local_recv_time += (end - begin);


		//only signal the odd ranks on other nodes that you have data - except master node
		MPI_Request request[NUM_CHANNELS];
		int i;
		for (i=0;i<NUM_CHANNELS;i++) {
			if 
			((interested_peers[i] != 0) && 
			(interested_peers[i] != 1)) {
				MPI_Isend(&alien_index, 1, MPI_INT, interested_peers[i], 2, *comm, &request[i]);
			}
			else {
				request[i] = MPI_REQUEST_NULL;
			}
		}
		MPI_Waitall(NUM_CHANNELS, request, MPI_STATUSES_IGNORE);

	}
}

/*
1. odd-numbered rank receives bit update from even-numbered rank on another node; if interested, it receives the data from this node;
2. odd-numbered rank also needs regularly to update the records of the even-numbered rank on localhost
*/
void *get_data(void *args) {
	MPI_Comm *comm = (MPI_Comm *) args;
	int rank;
	MPI_Comm_rank(*comm, &rank);
	MPI_Status status;
	int alien_index;
	//MPI_Request req;
	
	while (1) {
		//bit update from even-numbered rank
		MPI_Recv(&alien_index, 1, MPI_INT, MPI_ANY_SOURCE, 2, *comm, &status);
		int interesting_bit = has_interesting_bits(alien_index, status.MPI_SOURCE);
		if (interesting_bit != 0) {
			//request data from peer
			//MPI_Isend(&alien_index, 1, MPI_INT, status.MPI_SOURCE, 3, comm, &req); 
			MPI_Send(&alien_index, 1, MPI_INT, status.MPI_SOURCE, 3, *comm); 
			//data update from even-numbered rank on another node
			double begin = MPI_Wtime();
			MPI_Recv(&local_buffer[alien_index *  CHUNK_SIZE],  CHUNK_SIZE , MPI_CHAR, status.MPI_SOURCE, 4, *comm, MPI_STATUS_IGNORE);
			//printf("DEBUG: got data : %d from %d, bit %d\n", rank, status.MPI_SOURCE, alien_index);
			double end = MPI_Wtime();
			remote_recv_time += (end - begin);
			//send data to local process responsible for sending to other nodes
			MPI_Send(&alien_index, 1, MPI_INT, (rank - 1), 5, *comm);
			MPI_Send(&local_buffer[alien_index * CHUNK_SIZE], CHUNK_SIZE, MPI_CHAR, (rank - 1), 4, *comm);
		}
	}
}

/* even-number rank receives request from odd-number rank on another node and sends data
*/
void *send_data(void *args) {
	MPI_Comm *comm = (MPI_Comm *) args;
	MPI_Status status;
	int alien_index;
	//MPI_Request req;
	while (1) {
		MPI_Recv(&alien_index, 1, MPI_INT, MPI_ANY_SOURCE, 3, *comm, &status);
		double begin = MPI_Wtime();
		//MPI_Isend(&local_buffer[alien_index * CHUNK_SIZE], CHUNK_SIZE, MPI_CHAR, status.MPI_SOURCE, 4, comm, &req);
		MPI_Send(&local_buffer[alien_index * CHUNK_SIZE], CHUNK_SIZE, MPI_CHAR, status.MPI_SOURCE, 4, *comm);
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
	mpi_data_t mpi_data;
	int required = MPI_THREAD_MULTIPLE;
	MPI_Init_thread(&argc, &argv, required, &provided);
	printf("DEBUG: after init\n");
	MPI_Comm comm = MPI_COMM_WORLD;
	//try to read a hosts file which explicitly gives the mapping process-host
	FILE *fp;
	char name[MPI_MAX_PROCESSOR_NAME];
	if (NULL != (fp = fopen("hosts", "r"))) {
		MPI_Comm newcomm;
		int len;
		MPI_Get_processor_name(name, &len);
		int color = 1;
		int key = 0;
		char line[MPI_MAX_PROCESSOR_NAME];
		while (fgets(line,sizeof(line),fp)) {
			if (strncmp(line, name,len) == 0) {
				break;
			}
			else {
				key++;
			}
		}
		MPI_Comm_split(comm, color, key, &newcomm);
		comm = newcomm;
	}
	else {
		printf("Couldn't find hosts file\n");
	}
	MPI_Comm_rank(comm, &rank);
	mpi_data.rank = rank;
	//printf("DEBUG: rank %d, hostname %s\n", rank, name);
	MPI_Comm_size(comm, &size);
	mpi_data.size = size;
	pthread_t send_data_thread, get_local_data_thread, get_data_thread, terminator_thread;
	mpi_data.comm = comm;
	if (provided !=  required) {
		printf("Provided = %d,required = %d \n", provided, required);
		MPI_Abort(comm, -1);
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

	time1 = MPI_Wtime();

	map(rank, size);
	map2(rank, size);

	if (rank == 0) {
		master_send_shuffled(comm);
	}
	//even numbers send to (uneven, other nodes) or get (even, local node)
	if (rank == 0)
		pthread_create(&send_data_thread, NULL, send_data, &comm);
	else if (rank % 2 == 0) {
		pthread_create(&get_local_data_thread, NULL, even_get_local_data, &comm);
		pthread_create(&send_data_thread, NULL, send_data, &comm);

		pthread_mutex_lock(&mutex_ready);
		pthread_cond_wait(&signal_ready, &mutex_ready);
		pthread_mutex_unlock(&mutex_ready);
	}
	//odd numbers get data from other nodes and send local updates
	else {
		pthread_create(&get_data_thread, NULL, get_data, &comm);
	}

	double time2 = MPI_Wtime() - time1;
	double max_time;
	MPI_Barrier(comm);
	MPI_Reduce(&time2, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	if (rank == 0) {
		printf("TOTAL TIME: %lf\n", max_time);
	}
	MPI_Reduce(&local_recv_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	if (rank == 0) {
		printf("LOCAl RECEIVE TIME: %lf\n", max_time);
	}

	MPI_Reduce(&remote_recv_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	if (rank == 0) {
		printf("REMOTE RECEIVE TIME: %lf\n", max_time);
	}

	MPI_Reduce(&send_timer, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	if (rank == 0) {
		printf("TOTAL SEND TIME: %lf\n", max_time);
	}
	//printf("process %d to cancel threads\n", mpi_data.rank);
	MPI_Finalize();
}
