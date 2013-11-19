/*!
 * \file
 * Test for CBLAS wrappers
 * \include cblas_test.c
 */
#include "config.h"
#include "fupermod_cblas.h"
#include "fupermod/fupermod_conf.h"

#include <getopt.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

int main(int argc, char** argv) {
  char* conf_file = "./conf_file";
  int m = 1000;
  int n = 1000;
  int k = 1000;
  int opt;
  while ((opt = getopt(argc, argv, "hm:n:k:")) != -1)
    switch(opt) {
      case 'h':
        fprintf(stderr,
                "Test for mxm_cblas\n"
                "Usage: mxm_cblas_test [options]\n"
                "	-h	help\n"
                "	-m I	m size (default: %d)\n"
                "	-n I	n size (default: %d)\n"
                "	-k I	k size (default: %d)\n",
                m, n, k);
        return 0;
      case 'm':
        m = atoi(optarg);
        break;
      case 'n':
        n = atoi(optarg);
        break;
      case 'k':
        k = atoi(optarg);
        break;
    }
  
  int file_exist = 0;
  if(access(conf_file, R_OK) !=-1) {
    file_exist = 1;
  }
  if(file_exist==0){
    FILE* file = fopen(conf_file, "w");
    if(file == NULL) {
      perror("can't create conf_file");
      MPI_Finalize();
      return 0;
    }
    fupermod_print_conf(MPI_COMM_WORLD, 0, file, "cpu", "");
    fclose(file);
    fflush(file);
  }
  // configuration
  fupermod_process_conf conf = fupermod_get_conf(MPI_COMM_WORLD, conf_file);  
  fupermod_gemm* gemm = fupermod_gemm_alloc(&conf);

  fupermod_float* A = (fupermod_float*)malloc(sizeof(fupermod_float) * m * k);
  fupermod_float* B = (fupermod_float*)malloc(sizeof(fupermod_float) * k * n);
  fupermod_float* C = (fupermod_float*)malloc(sizeof(fupermod_float) * m * n);
  int i;
  for (i = 0; i < m * k; i++)
    A[i] = 1.1;
  for (i = 0; i < k * n; i++)
    B[i] = 2.5;
  for (i = 0; i < m * n; i++)
    C[i] = 0.0;

  struct timeval start, end;
  gettimeofday(&start, NULL);
  fupermod_gemm_execute(gemm, CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1,
                A, k, B, n, 0, C, n);
  gettimeofday(&end, NULL);
  double time = ((end.tv_sec + end.tv_usec / 1000000.) -
                 (start.tv_sec + start.tv_usec / 1000000.));
  printf("#m\tn\tk\ttime\n");
  printf("%d\t%d\t%d\t%le\n", m, n, k, time);

  free(A);
  free(B);
  free(C);

  fupermod_gemm_free(gemm);
  return 0;
}
