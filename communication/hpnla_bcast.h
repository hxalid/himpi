#ifndef HPNLA_BCAST_H_
#define HPNLA_BCAST_H_

#include <mpi.h>

typedef enum hpnla_bcast_algo {
  original,   //0
  lin,
  binary,
  flat,
  binomial,
  scatter_lr_allgather,
  scatter_rd_allgather,
  binomial_mpich,
  pipelined_linear // 8
} hpnla_bcast_algo;

void hpnla_bcast(void *buffer, int count, MPI_Datatype datatype,
    int root, MPI_Comm comm, hpnla_bcast_algo algorithm);

void bcast_scatter_lr_allgather(void * buff, int count, MPI_Datatype data_type,
		int root, MPI_Comm comm);


#endif /*HPNLA_BCAST_H_*/
