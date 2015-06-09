#include "hpnla_cblas.h"
#include <mpi.h>
double two_dot_five(hpnla_gemm* gemm,
                    size_t m,
                    size_t k,
                    size_t n,
                    size_t Block_size,
                    size_t group,
                    size_t key,
                    size_t size_row,
                    size_t size_col,
                    size_t NB_groups,
                    MPI_Comm world,
                    int bcast_algorithm,
                    int distribution
                    );

