#include "hpnla_bcast.h"
#include "tools/hpnla_debug.h"
#include <stdio.h>


#define Bcast_TAG 904920477 
void hpnla_bcast(void *buffer, int count, MPI_Datatype datatype,
    int root, MPI_Comm comm, hpnla_bcast_algo algorithm){
  int myrank;
  int size;
  MPI_Status status;
  MPI_Comm_rank ( comm, &myrank );
  MPI_Comm_size ( comm, &size );
  if (size ==1 ) return;
  int pos = (myrank - root+size) %size;
  int i=0;
  myrank +=size;
  switch (algorithm){
    case lin:
      if(pos == 0){
        debug_print_rank(2, (size_t)myrank, "Sent Root : %d myrank : %d    \
            sent to rank: %d mypos: %d sent to pos : %d\n",
            root, myrank, (myrank+1)%size, pos, (pos+1)%size);
        MPI_Send(buffer, count, datatype, (myrank+1)%size,
            Bcast_TAG, comm);
      }else{
        debug_print_rank(2, (size_t)myrank, "Recv Root : %d myrank : %d    \
            recv from rank: %d mypos: %d recv from pos : %d\n",
            root, myrank, (myrank-1)%size, pos, (pos-1)%size);
        MPI_Recv(buffer, count, datatype, (myrank-1+size)%size,
            Bcast_TAG, comm, &status);
        if((pos+1) < size)
        {
          debug_print_rank(2, (size_t)myrank, "Sent Root : %d myrank : %d    \
              sent to rank: %d mypos: %d sent to pos : %d\n",
              root, myrank, (myrank+1)%size, pos, (pos+1)%size);
          MPI_Send(buffer, count, datatype, (myrank+1)%size,
              Bcast_TAG, comm);
        }
      }
      break;
    case binary:
      if(pos != 0){
        int rec_pos = (root + (int)((pos-1)/2)+size) %size;
        debug_print_rank(2, (size_t)myrank, "Recv Root : %d myrank : %d    \
            recv from rank: %d mypos: %d recv from pos : %d\n",
            root, myrank, rec_pos, pos, (int)((pos-1)/2));
        MPI_Recv(buffer, count, datatype, rec_pos, Bcast_TAG, comm, &status);
      }
      if(pos * 2 + 1 < size)
      {
        int send_pos = (root + pos * 2 + 1)%size;
        debug_print_rank(2, (size_t)myrank, "Sent Root : %d myrank : %d    \
            sent to rank: %d mypos: %d sent to pos : %d\n",
            root, myrank, send_pos, pos, (pos * 2 + 1)%size);
        MPI_Send(buffer, count, datatype, send_pos,
            Bcast_TAG, comm);
      }
      if(pos * 2 + 2 < size)
      {
        int send_pos = (root + pos * 2 + 2 +size)%size;
        debug_print_rank(2, (size_t)myrank, "Sent Root : %d myrank : %d    \
            sent to rank: %d mypos: %d sent to pos : %d\n",
            root, myrank, send_pos, pos, (pos * 2 + 2)%size);
        MPI_Send(buffer, count, datatype, send_pos, Bcast_TAG, comm);
      }
      break;
    case binomial:
      for(i=1; i < size; i = i*2){
        if(pos < i){
          if(pos+i <size)
            debug_print_rank(2, (size_t)myrank, "Sent Root : %d myrank : %d    \
                sent to rank: %d mypos: %d sent to pos : %d\n",
                root, myrank, (myrank + i)%size, pos, pos+i);
          MPI_Send(buffer, count, datatype, (myrank + i)%size,
              Bcast_TAG, comm);
        }else if(pos < i*2){
          debug_print_rank(2, (size_t)myrank, "Recv Root : %d myrank : %d    \
              recv from rank: %d mypos: %d recv from pos : %d\n",
              root, myrank, (myrank - i)%size, pos,
              ((pos-i+size)%size));
          MPI_Recv(buffer, count, datatype, (myrank - i)%size,
              Bcast_TAG, comm, &status);
        }

      }
      break;
    case flat:
      if(myrank % size == root){
        int i = 0;
        for (i = 0; i < size; i++)
          if(i!=root)
          {
            debug_print_rank(2, (size_t)myrank, "Sent Root : %d myrank : %d    \
                sent to rank: %d mypos: %d sent to pos : %d\n",
                root, myrank,  i, pos, 0);
            MPI_Send(buffer, count, datatype, i,
                Bcast_TAG, comm);
          }
      }else{
        debug_print_rank(2, (size_t)myrank, "Recv Root : %d myrank : %d    \
            recv from rank: %d mypos: %d recv from pos : %d\n",
            root, myrank, root, pos,
            0);
        MPI_Recv(buffer, count, datatype, root,
            Bcast_TAG, comm, &status);
      }
      break;
    case original:
      MPI_Bcast(buffer, count, datatype, root, comm);
      break;
    case scatter_lr_allgather:
      bcast_scatter_lr_allgather(buffer, count, datatype, root, comm);
      break;
    default:
      fprintf(stderr, "algorithm unknow we used the MPI one\n");
      MPI_Bcast(buffer, count, datatype, root, comm);
      break;
  }
}

/*
 * Originally implemented inside MPICH. Here it has been modified.
 *
 * A ring algorithm is used for the allgather, which takes p-1 steps.
 * This may perform better than recursive doubling for long messages and
 * medium-sized non-power-of-two messages.
 * Total Cost = (lgp+p-1).alpha + 2.n.((p-1)/p).beta
 */
void bcast_scatter_lr_allgather(void * buff, int count, MPI_Datatype data_type,
		int root, MPI_Comm comm) {

	MPI_Aint extent;
		MPI_Status status;
		int i, j, k, src, dst, rank, num_procs;
		int mask, relative_rank, curr_size, recv_size, send_size, nbytes;
		int scatter_size, left, right, next_src, *recv_counts, *disps;
		int tag = 1, success = 0, failure = 1;

		/*	MPI_Comm comm = _comm;
		 MPI_Comm_dup(_comm, &comm); */

		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &num_procs);

		int data_type_size;
		MPI_Type_size(data_type, &data_type_size);

		nbytes = data_type_size*count;

		scatter_size = (nbytes + num_procs - 1) / num_procs; 	// ceiling division
		curr_size = (rank == root) ? nbytes : 0; 				// root starts with all the data
		relative_rank = (rank >= root) ? rank - root : rank - root + num_procs;

		mask = 0x1;     // start from 1
		while (mask < num_procs) {
			if (relative_rank & mask) {
				src = rank - mask;
				if (src < 0) {
					src += num_procs;
				}

				recv_size = nbytes - relative_rank * scatter_size;

				//  recv_size is larger than what might actually be sent by the
				//  sender. We don't need compute the exact value because MPI
				//  allows you to post a larger recv.
				if (recv_size <= 0) {
					curr_size = 0; // this process doesn't receive any data because of uneven division
				} else {
					MPI_Recv(buff + relative_rank * scatter_size, recv_size,
							MPI_BYTE, src, tag, comm, &status);
					MPI_Get_count(&status, MPI_BYTE, &curr_size);
				}

				break;
			}
			mask <<= 1;  // mask = mask*2
		}

		// This process is responsible for all processes that have bits
		// set from the LSB upto (but not including) mask.  Because of
		// the "not including", we start by shifting mask back down
		// one.

		mask >>= 1;
		while (mask > 0) {
			if (relative_rank + mask < num_procs) {
				send_size = curr_size - scatter_size * mask;
				// mask is also the size of this process's subtree

				if (send_size > 0) {
					dst = rank + mask;
					if (dst >= num_procs)
						dst -= num_procs;

					MPI_Send(buff + scatter_size * (relative_rank + mask),
							send_size, MPI_BYTE, dst, tag, comm);

					curr_size -= send_size;
				}
			}
			mask >>= 1;  // mask = mask/2
		}

		// done scatter now do allgather
		recv_counts = (int *) malloc(sizeof(int) * num_procs);
		if (!recv_counts) {
			fprintf(stderr, "bcast-scatter-LR-allgather: cannot allocate memory for recv_counts\n");
			MPI_Finalize();
			exit(failure);
		}

		disps = (int *) malloc(sizeof(int) * num_procs);
		if (!disps) {
			fprintf(stderr, "bcast-scatter-LR-allgather: cannot allocate memory for disps\n");
			MPI_Finalize();
			exit(failure);
		}

		for (i = 0; i < num_procs; i++) {
			recv_counts[i] = nbytes - i * scatter_size;
			if (recv_counts[i] > scatter_size)
				recv_counts[i] = scatter_size;
			if (recv_counts[i] < 0)
				recv_counts[i] = 0;
		}

		disps[0] = 0;
		for (i = 1; i < num_procs; i++)
			disps[i] = disps[i - 1] + recv_counts[i - 1];

		left = (num_procs + rank - 1) % num_procs;
		right = (rank + 1) % num_procs;

		src = rank;
		next_src = left;

		/*
		 * Now lets do all gather on ring
		 */
		for (i = 1; i < num_procs; i++) {
			MPI_Sendrecv(buff + disps[(src - root + num_procs) % num_procs],
					recv_counts[(src - root + num_procs) % num_procs], MPI_BYTE,
					right, tag,
					buff + disps[(next_src - root + num_procs) % num_procs],
					recv_counts[(next_src - root + num_procs) % num_procs],
					MPI_BYTE, left, tag, comm, &status);

			src = next_src;
			next_src = (num_procs + next_src - 1) % num_procs;
		}

		free(recv_counts);
		free(disps);
}



