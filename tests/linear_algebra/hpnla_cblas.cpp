#include "config.h"
#include "hpnla_cblas.h"


#ifdef HAVE_LIBMKL_CORE
#include <mkl_blas.h>
#include <mkl_service.h>
#endif

#if defined HAVE_LIBACML || defined HAVE_LIBACML_MP
#include <omp.h>
#include <acml.h>
#endif

#ifdef HAVE_LIBESSLBG
#include <essl.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern "C"
hpnla_gemm* hpnla_gemm_alloc(hpnla_process_conf* conf) {
  hpnla_gemm* gemm = NULL;
  if (conf) {
    gemm = (hpnla_gemm*)malloc(sizeof(hpnla_gemm));
    gemm->conf = conf;
    gemm->time = 0;
    // process binding, also suitable for a process with multiple threads
    if (strcmp(conf->bind, "na")) {
      hpnla_bind_process(conf->bind);
    }
    // parse suboptions
    char* subopts = strdup(conf->subopts);
    char* subopts_0 = subopts;
    if (!strcmp(conf->device_type, "cpu")) {
      char* tokens[] = {"t", NULL};
      char* value;
      while (*subopts != '\0') {
        switch (getsubopt(&subopts, tokens, &value)) {
          case 0:
            int threads = atoi(value);
            #ifdef HAVE_LIBMKL_CORE
            mkl_set_num_threads(threads);
            #elif defined HAVE_LIBACML_MP
            if (strcmp(conf->bind, "na")){
              omp_set_num_threads(threads);
            }
            #endif
            break;
        }
      }
      gemm->params = NULL;
    } 
    free(subopts_0);
  }
  return gemm;
}

char CBLAS_TRANSPOSE_CHAR[] = {'N', 'T', 'C'};
char* cblas_transpose(CBLAS_TRANSPOSE Trans) {
  switch(Trans) {
    case 111: return &CBLAS_TRANSPOSE_CHAR[0];
    case 112: return &CBLAS_TRANSPOSE_CHAR[1];
    case 113: return &CBLAS_TRANSPOSE_CHAR[2];
  }
  return NULL;
}

extern	"C"
void hpnla_gemm_execute(struct hpnla_gemm* gemm,
		const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
		const int M, const int N, const int K,
		const hpnla_float alpha, const hpnla_float *A, const int lda,
		const hpnla_float *B, const int ldb,  const hpnla_float beta,
		hpnla_float *C, const int ldc) 
{
  if (!gemm || !strcmp(gemm->conf->device_type, "cpu")) { // !gemm: in case conf = NULL and then gemm = NULL
      #ifdef WITH_BLAS_GSL
      cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      #elif defined HAVE_LIBGOTO2
      cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      #elif defined HAVE_LIBMKL_CORE
      dgemm(cblas_transpose(TransB), cblas_transpose(TransA), &N, &M, &K, &alpha, B, &ldb, A, &lda, &beta, C, &ldc);
      #elif defined HAVE_LIBACML || defined HAVE_LIBACML_MP
      //note: BLAS is colum-major, so have to swap matrices A and B
      dgemm(*cblas_transpose(TransB), *cblas_transpose(TransA), N, M, K, alpha, (hpnla_float *)B, ldb, (hpnla_float *)A, lda, beta, C, ldc);
	    #elif defined HAVE_LIBESSLBG
      //note: BLAS is colum-major, so have to swap matrices A and B
      dgemm((char*)TransB, (char*)TransA, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
      #else
      dgemm_((char)TransB, (char)TransA, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
      #endif
  } 
}

extern "C"
void hpnla_gemm_free(hpnla_gemm* gemm) {
  if (gemm) {
    if (gemm->params) {
      free(gemm->params);
    }
    free(gemm);
  }
}
