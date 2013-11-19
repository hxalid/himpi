#include "hdnla_cblas.h"
#include <mpi.h>
double cannon(hdnla_gemm  * gemm,
                   double *a, double *b, double *c,
                   size_t local_matrix_dim,int proc_assign,int cart_coords[],MPI_Comm comm_cart);


