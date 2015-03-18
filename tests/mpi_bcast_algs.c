#include "../tests/mpi_bcast_algs.h"



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

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_procs);

    int data_type_size;
    MPI_Type_size(data_type, &data_type_size);

    nbytes = data_type_size*count;

    scatter_size = (nbytes + num_procs - 1) / num_procs; // ceiling division
    curr_size = (rank == root) ? nbytes : 0; // root starts with all the data
    relative_rank = (rank >= root) ? rank - root : rank - root + num_procs;

    mask = 0x1; // start from 1
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
        mask <<= 1; // mask = mask*2
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
        mask >>= 1; // mask = mask/2
    }

    // done scatter now do allgather
    recv_counts = (int *) malloc(sizeof (int) * num_procs);
    if (!recv_counts) {
        fprintf(stderr, "bcast-scatter-LR-allgather: cannot allocate memory for recv_counts\n");
        MPI_Finalize();
        exit(failure);
    }

    disps = (int *) malloc(sizeof (int) * num_procs);
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

/*****************************************************************************

 * Function: bcast_scatter_rdb_allgather

 * Return: int

 * Inputs:
    buff: send input buffer
    count: number of elements to send
    data_type: data type of elements being sent
    root: source of data
    comm: communicator

 * Descrp: broadcasts using a scatter followed by rdb allgather.

 * Auther: MPICH / modified by Ahmad Faraj

 ****************************************************************************/

int bcast_scatter_rdb_allgather(void * buff, int count, MPI_Datatype data_type,
        int root, MPI_Comm comm) {
    MPI_Aint extent;
    MPI_Status status;

    int i, j, k, src, dst, rank, num_procs, send_offset, recv_offset;
    int mask, relative_rank, curr_size, recv_size, send_size, nbytes;
    int scatter_size, tree_root, relative_dst, dst_tree_root;
    int my_tree_root, offset, tmp_mask, num_procs_completed;
    int tag = 1, success = 0, failure = 1;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_procs);
    MPI_Type_extent(data_type, &extent);

    nbytes = extent * count;
    scatter_size = (nbytes + num_procs - 1) / num_procs; // ceiling division 
    curr_size = (rank == root) ? nbytes : 0; // root starts with all the data
    relative_rank = (rank >= root) ? rank - root : rank - root + num_procs;

    mask = 0x1;
    while (mask < num_procs) {
        if (relative_rank & mask) {
            src = rank - mask;
            if (src < 0) src += num_procs;
            recv_size = nbytes - relative_rank * scatter_size;
            //  recv_size is larger than what might actually be sent by the
            //  sender. We don't need compute the exact value because MPI
            //  allows you to post a larger recv.
            if (recv_size <= 0)
                curr_size = 0; // this process doesn't receive any data
                // because of uneven division 
            else {
                MPI_Recv(buff + relative_rank * scatter_size, recv_size,
                        MPI_BYTE, src, tag, comm, &status);
                MPI_Get_count(&status, MPI_BYTE, &curr_size);
            }
            break;
        }
        mask <<= 1;
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
                if (dst >= num_procs) dst -= num_procs;
                MPI_Send(buff + scatter_size * (relative_rank + mask),
                        send_size, MPI_BYTE, dst, tag, comm);

                curr_size -= send_size;
            }
        }
        mask >>= 1;
    }

    // done scatter now do allgather

    mask = 0x1;
    i = 0;
    while (mask < num_procs) {
        relative_dst = relative_rank ^ mask;

        dst = (relative_dst + root) % num_procs;

        /* find offset into send and recv buffers.
           zero out the least significant "i" bits of relative_rank and
           relative_dst to find root of src and dst
           subtrees. Use ranks of roots as index to send from
           and recv into  buffer */

        dst_tree_root = relative_dst >> i;
        dst_tree_root <<= i;

        my_tree_root = relative_rank >> i;
        my_tree_root <<= i;

        send_offset = my_tree_root * scatter_size;
        recv_offset = dst_tree_root * scatter_size;

        if (relative_dst < num_procs) {
            MPI_Sendrecv(buff + send_offset, curr_size, MPI_BYTE, dst, tag,
                    buff + recv_offset, scatter_size * mask, MPI_BYTE, dst,
                    tag, comm, &status);
            MPI_Get_count(&status, MPI_BYTE, &recv_size);
            curr_size += recv_size;
        }

        /* if some processes in this process's subtree in this step
           did not have any destination process to communicate with
           because of non-power-of-two, we need to send them the
           data that they would normally have received from those
           processes. That is, the haves in this subtree must send to
           the havenots. We use a logarithmic recursive-halfing algorithm
           for this. */

        if (dst_tree_root + mask > num_procs) {
            num_procs_completed = num_procs - my_tree_root - mask;
            /* num_procs_completed is the number of processes in this
               subtree that have all the data. Send data to others
               in a tree fashion. First find root of current tree
               that is being divided into two. k is the number of
               least-significant bits in this process's rank that
               must be zeroed out to find the rank of the root */
            j = mask;
            k = 0;
            while (j) {
                j >>= 1;
                k++;
            }
            k--;

            offset = scatter_size * (my_tree_root + mask);
            tmp_mask = mask >> 1;

            while (tmp_mask) {
                relative_dst = relative_rank ^ tmp_mask;
                dst = (relative_dst + root) % num_procs;

                tree_root = relative_rank >> k;
                tree_root <<= k;

                /* send only if this proc has data and destination
                   doesn't have data. */

                if ((relative_dst > relative_rank)
                        && (relative_rank < tree_root + num_procs_completed)
                        && (relative_dst >= tree_root + num_procs_completed)) {
                    MPI_Send(buff + offset, recv_size, MPI_BYTE, dst, tag, comm);

                    /* recv_size was set in the previous
                       receive. that's the amount of data to be
                       sent now. */
                }/* recv only if this proc. doesn't have data and sender
		 has data */
                else if ((relative_dst < relative_rank)
                        && (relative_dst < tree_root + num_procs_completed)
                        && (relative_rank >= tree_root + num_procs_completed)) {

                    MPI_Recv(buff + offset, scatter_size * num_procs_completed,
                            MPI_BYTE, dst, tag, comm, &status);

                    /* num_procs_completed is also equal to the no. of processes
                       whose data we don't have */
                    MPI_Get_count(&status, MPI_BYTE, &recv_size);
                    curr_size += recv_size;
                }
                tmp_mask >>= 1;
                k--;
            }
        }
        mask <<= 1;
        i++;
    }

    return success;
}

/*
 * pipelined linear tree
 */
int bcast_pipelined_linear(void *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
    int tag = 5000;
    MPI_Status status;
    MPI_Request req;
    int rank, phase, pipe_length;
    int i, j;
    int last_count;
    int send_size;
    int total_node;

    int *send_to;
    int *recv_from;
    int *sequence;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_node);

    // TOPOLOGY UNAWARE                                                                                                           
    send_to = (int *) malloc(total_node * sizeof (int));
    recv_from = (int *) malloc(total_node * sizeof (int));
    sequence = (int *) malloc(total_node * sizeof (int));

    for (i = 0; i < total_node; i++) {
        send_to[i] = i + 1;
        recv_from[i] = i - 1;
        sequence[i] = i;
    }
    recv_from[0] = -1;
    send_to[total_node - 1] = 0;

    // when a message is smaller than a block size 
    if (count <= bcast_linear_segment_size) {
        for (i = 0; i < total_node - 1; i++) {
            if ((send_to[rank] != -1) && (sequence[rank] == i))
                MPI_Send(buf, count, datatype, send_to[rank], tag, MPI_COMM_WORLD);

            if ((recv_from[rank] != -1) && (sequence[rank] == i + 1))
                MPI_Recv(buf, count, datatype, recv_from[rank], tag, MPI_COMM_WORLD, &status);
        }

        free(send_to);
        free(recv_from);
        free(sequence);
        
       // if (!rank)
           printf("rank=%d, step 0\n", rank);
        
        return 0;
    } 
    
    pipe_length = (count - 1) / bcast_linear_segment_size + 1;
    last_count = count % bcast_linear_segment_size;
    phase = pipe_length + total_node - 2;

    // last message has the same size as all the rest
    if (last_count == 0) {
        for (i = 0; i < phase; i++) {
            if ((send_to[rank] != -1) && (i >= sequence[rank]) && (i < (pipe_length + sequence[rank])))
                MPI_Send(buf + (bcast_linear_segment_size * (i - sequence[rank])), bcast_linear_segment_size, datatype, send_to[rank], tag, MPI_COMM_WORLD);

            if ((recv_from[rank] != -1) && (i >= sequence[rank] - 1) && (i < pipe_length + sequence[rank] - 1))
                MPI_Recv(buf + (bcast_linear_segment_size * (i - sequence[rank] + 1)), bcast_linear_segment_size, datatype, recv_from[rank], tag, MPI_COMM_WORLD, &status);
        }
    } else {
        send_size = bcast_linear_segment_size;
        for (i = 0; i < phase; i++) {
            if ((send_to[rank] != -1) && (i >= sequence[rank]) && (i < (pipe_length + sequence[rank]))) {
                if (i == (pipe_length + sequence[rank] - 1))
                    send_size = last_count;
                MPI_Send(buf + (bcast_linear_segment_size * (i - sequence[rank])), send_size, datatype, send_to[rank], tag, MPI_COMM_WORLD);
            }
            if ((recv_from[rank] != -1) && (i >= sequence[rank] - 1) && (i < pipe_length + sequence[rank] - 1)) {
                if (i == (pipe_length + sequence[rank] - 2))
                    send_size = last_count;
                MPI_Recv(buf + (bcast_linear_segment_size * (i - sequence[rank] + 1)), send_size, datatype, recv_from[rank], tag, MPI_COMM_WORLD, &status);
            }
        }
    }

    free(send_to);
    free(recv_from);
    free(sequence);

    return 0;
}

/*****************************************************************************

 * Function: bcast_binomial_tree

 * Return: int

 * Inputs:
    buff: send input buffer
    count: number of elements to send
    data_type: data type of elements being sent
    root: source of data
    comm: communicator

 * Descrp: broadcasts using a binomial tree.

 * Auther: MPICH / modified by Ahmad Faraj

 ****************************************************************************/

int bcast_binomial_tree(void * buff, int count, MPI_Datatype data_type,
        int root, MPI_Comm comm) {
    MPI_Status status;
    int i, src, dst, rank, num_procs, mask, relative_rank;
    int tag = 1, success = 0;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_procs);

    relative_rank = (rank >= root) ? rank - root : rank - root + num_procs;

    //  printf("rank=%d, relative_rank = %d\n", rank, relative_rank);

    mask = 0x1;
    while (mask < num_procs) {
        if (relative_rank & mask) {
            src = rank - mask;
            if (src < 0) src += num_procs;
            //   printf("rank=%d, relative_rank=%d, mask=%d, src=%d\n", rank, relative_rank, mask, src);
            MPI_Recv(buff, count, data_type, src, tag, comm, MPI_STATUS_IGNORE);
            break;
        }
        mask <<= 1; // mask *= 2
    }

    mask >>= 1;
    while (mask > 0) {
        if (relative_rank + mask < num_procs) {
            dst = rank + mask;
            if (dst >= num_procs) dst -= num_procs;
            MPI_Send(buff, count, data_type, dst, tag, comm);
            //   printf("rank=%d, relative_rank=%d, mask=%d, dst=%d\n", rank, relative_rank, mask, dst);
        }
        mask >>= 1; // mask /= 2
    }

    return success;
}

void hpnla_bcast(void *buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm, hpnla_bcast_algo algorithm) {
    int myrank;
    int size;
    MPI_Status status;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &size);
    if (size == 1) return;
    int pos = (myrank - root + size) % size;
    int i = 0;
    myrank += size;
    switch (algorithm) {
        case lin:
            if (pos == 0) {
                MPI_Send(buffer, count, datatype, (myrank + 1) % size,
                        Bcast_TAG, comm);
            } else {
                MPI_Recv(buffer, count, datatype, (myrank - 1 + size) % size,
                        Bcast_TAG, comm, &status);
                if ((pos + 1) < size) {
                    MPI_Send(buffer, count, datatype, (myrank + 1) % size,
                            Bcast_TAG, comm);
                }
            }
            break;
        case binary:
            if (pos != 0) {
                int rec_pos = (root + (int) ((pos - 1) / 2) + size) % size;
                MPI_Recv(buffer, count, datatype, rec_pos, Bcast_TAG, comm, &status);
            }
            if (pos * 2 + 1 < size) {
                int send_pos = (root + pos * 2 + 1) % size;
                MPI_Send(buffer, count, datatype, send_pos,
                        Bcast_TAG, comm);
            }
            if (pos * 2 + 2 < size) {
                int send_pos = (root + pos * 2 + 2 + size) % size;
                MPI_Send(buffer, count, datatype, send_pos, Bcast_TAG, comm);
            }
            break;
        case binomial:
            for (i = 1; i < size; i = i * 2) {
                if (pos < i) {
                    if (pos + i < size)
                        MPI_Send(buffer, count, datatype, (myrank + i) % size,
                            Bcast_TAG, comm);
                } else if (pos < i * 2) {
                    MPI_Recv(buffer, count, datatype, (myrank - i) % size,
                            Bcast_TAG, comm, &status);
                }

            }
            break;
        case flat:
            if (myrank % size == root) {
                int i = 0;
                for (i = 0; i < size; i++)
                    if (i != root) {
                        MPI_Send(buffer, count, datatype, i,
                                Bcast_TAG, comm);
                    }
            } else {
                MPI_Recv(buffer, count, datatype, root,
                        Bcast_TAG, comm, &status);
            }
            break;
        case original: //4
            MPI_Bcast(buffer, count, datatype, root, comm);
            break;
        case scatter_lr_allgather:
            bcast_scatter_lr_allgather(buffer, count, datatype, root, comm);
            break;
        case scatter_rd_allgather:
            bcast_scatter_rdb_allgather(buffer, count, datatype, root, comm);
            break;
        case binomial_mpich:
            bcast_binomial_tree(buffer, count, datatype, root, comm);
            break;
        case pipelined_linear:
            bcast_pipelined_linear(buffer, count, datatype, root, comm);
            break;
        default:
            fprintf(stderr, "Unknown algorithm we used the MPI one\n");
            MPI_Bcast(buffer, count, datatype, root, comm);
            break;
    }
}


