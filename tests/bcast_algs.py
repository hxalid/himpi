#!/usr/bin/python

root = 0
size = 32
rec_pos = 0
send_pos = 0

for myrank in range(0, 32):
    pos = (myrank - root + size) % size

    if (pos != 0):
        rec_pos = (root + (pos - 1) / 2 + size) % size
        print("REC: %s <- %s, pos = %s" % (myrank, rec_pos, pos))
        #MPI_Recv(buffer, count, datatype, rec_pos, Bcast_TAG, comm, &status);
    

    if (pos * 2 + 1 < size):
        send_pos = (root + pos * 2 + 1) % size
        #MPI_Send(buffer, count, datatype, send_pos, Bcast_TAG, comm)
        print("SEND: %s -> %s, pos = %s" % (myrank, send_pos, pos))
    
    
    if (pos * 2 + 2 < size):
        send_pos = (root + pos * 2 + 2 + size) % size
        #MPI_Send(buffer, count, datatype, send_pos, Bcast_TAG, comm)
        print("SEND1: %s -> %s, pos = %s" % (myrank, send_pos, pos))

#    print("myrank = %s, pos = %s, rec_pos = %s, send_post = %s" % (myrank, pos, rec_pos, send_pos))
    



