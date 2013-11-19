/*!
 * \file
 * Block Matrix Multiplication example
 *
 * Authors: Quintin Jean-Noël
 */

#include "config.h"
#include "hpnla_cblas.h"
#include "hpnla_memory.h"
#include "Matrix_init.h"

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#ifdef HDNLA_SMPI
#include <smpi.h>
#else
#define SMPI_SHARED_MALLOC malloc
#define SMPI_SHARED_FREE free
#endif



int main(int argc, char ** argv){
  double *a  , *b  , *c  ; //matrices
  double *B_a, *B_b, *B_c; //matrix blocks
  double *check;           //matrices for checking the result
  double alpha=1, beta=1;  //C := alpha * a * b + beta * c

  size_t    m   = 1000, n   = 1000, k   = 1000; // matrix sizes
  size_t    B_m = 100 , B_n = 100 , B_k = 100 ; // block sizes
  size_t    lda = k   , ldb = n   , ldc = n   ; // matrix line size
  size_t    x   = 0   , y   = 0   , z   = 0   ; // iterators
  /* x index on M
     y index on N
     Z index on K */



  int opt;
  optind = 1;

  //get the parameter from command line
  while ((opt = getopt(argc, argv, "hM:N:K:m:n:k:")) != -1)
    switch(opt) {
      case 'h':
        fprintf(stderr,
                "Usage: mxm_cblas_test [options]\n"
                "	-M I	M size (default: %zu)\n"
                "	-N I	N size (default: %zu)\n"
                "	-K I	K size (default: %zu)\n"
                "	-m I	B_M block size (default: %zu)\n"
                "	-n I	B_N block size (default: %zu)\n"
                "	-k I	B_K block size (default: %zu)\n"
                "	-h	help\n",
                m, n, k, B_m, B_n, B_k);
        return 0;
      case 'M':
        m = atoi(optarg);
        break;
      case 'N':
        n   = atoi(optarg);
        ldb = n;
        ldc = n;
        break;
      case 'K':
        k   = atoi(optarg);
        lda = k;
        break;
      case 'm':
        B_m = atoi(optarg);
        break;
      case 'n':
        B_n = atoi(optarg);
        break;
      case 'k':
        B_k = atoi(optarg);
        break;
    }

  // Defined the device if we use the GPU
  //TODO explain parameters
  hpnla_gemm  * gemm = hpnla_gemm_alloc(NULL); // Allocate the gemm workspace


  //TODO explain the parameter HDNLA_MEM_DEFAULT
  a     =  (double *) SMPI_SHARED_MALLOC(sizeof(double) * m * k);
  if ( a == 0 ){
    perror("Error allocation Matrix A");
    exit(-1);
  }

  //TODO explain the parameter HDNLA_MEM_DEFAULT
  b     = (double *) SMPI_SHARED_MALLOC(sizeof(double) * k * n);
  if ( b == 0 ){
    perror("Error allocation Matrix B");
    exit(-1);
  }

  //TODO explain the parameter HDNLA_MEM_DEFAULT
  c     = (double *) SMPI_SHARED_MALLOC(sizeof(double) * m * n);
  if ( c == 0 ){
    perror("Error allocation Matrix C");
    exit(-1);
  }

  //TODO explain the parameter HDNLA_MEM_DEFAULT
  check = (double *) SMPI_SHARED_MALLOC(sizeof(double) * m * n);
  if ( check == 0){
    perror("Error allocation Matrix Check");
    exit(-1);
  }

  // intialisation of the matrices
#if 0
  for( x=0; x<m; x++){
    for( z=0; z<k; z++){
      a[x*lda+z] = (double)(x+z);
    }
  }
  for( z=0; z<k; z++){
    for( y=0; y<n; y++){
      b[z*ldb+y] = (double)(y);
    }
  }
  for( x=0; x<m; x++){
    for( y=0; y<n; y++){
      c[x*ldc+y] = (double)0;
      check[x*ldc+y] = (double)0;
    }
  }
#else
  matrices_initialisation(&a, &b, &c, m, k, k, n, 0, 0);
#endif

  struct timeval start, end; //time mesure
  double time;
  gettimeofday(&start, NULL);

  // computation of the gemm
  for( x=0; x < m; x+=B_m){
    for( y=0; y < n; y+=B_n){
      for( z=0; z < k; z+=B_k){
        //first element of a block is defined by:
        // (x, z) or (z, y) or (x, y)
        B_a = a + x * lda + z ;
        B_b = b + z * ldb + y ;
        B_c = c + x * ldc + y ;

        hpnla_gemm_execute(gemm, //user parameters: here it's equal to NULL
                      CblasRowMajor, CblasNoTrans, CblasNoTrans,
                      B_m, B_n, B_k, alpha, B_a, lda, B_b, ldb,
                      beta, B_c, ldc);
      }
    }
  }

  gettimeofday(&end, NULL);
  time = ((end.tv_sec + end.tv_usec / 1000000.)
          - (start.tv_sec + start.tv_usec / 1000000.));
  printf("#m\tn\tk\tB_m\tB_n\tB_k\ttime\n");
  printf("%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%le\n", m, n, k, B_m, B_n, B_k, time);

  /* inital code for dgemm */

  gettimeofday(&start, NULL);

  hpnla_gemm_execute(gemm, //user parameters: here it's equal to NULL
                CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, alpha, a, lda, b, ldb,
                beta, check, ldc);

  gettimeofday(&end, NULL);
  time = ((end.tv_sec + end.tv_usec / 1000000.)
          - (start.tv_sec + start.tv_usec / 1000000.));
  printf("#m\tn\tk\tB_m\tB_n\tB_k\ttime\n");
  printf("%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%le\n", m, n, k, m, n, k, time);
#if 0
  /*Display for checking */
  bool result_good=true;
  for( x=0; x<m; x++){
    for( y=0; y<n; y++){
      /* WARNING this could be lead to some errors ( precision with double )*/
      if (c[x*ldc + y] != check[x*ldc + y]){
        result_good = false;
        printf("%lf\t%lf\n",c[x*ldc+y],check[x*ldc+y]);
      }
    }
  }
  if (result_good)
    printf("result check ok\n");
  else
    printf("result check not ok\n"
           "WARNING the test could be lead to some "
           "errors ( precision with double )\n");
#else
  check_result(c, a, b, m, n, k, k, 0, 0, 1, 1);
#endif
  // close properly the pragram
  SMPI_SHARED_FREE(a);
  SMPI_SHARED_FREE(b);
  SMPI_SHARED_FREE(c);
  SMPI_SHARED_FREE(check);

  hpnla_gemm_free(gemm);
  return 0;
}
