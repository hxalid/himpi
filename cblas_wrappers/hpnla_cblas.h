#ifndef HDNLA_CBLAS_H_
#define HDNLA_CBLAS_H_

#include "config.h"
#include "tools/hpnla_conf.h"
#include "tools/hpnla_bind.h"
#include <gsl/gsl_cblas.h>

/*!
 * \defgroup hpnla_cblas Wrappers for BLAS libraries
 * Main data structures are BLAS operation workspaces, for example \ref hpnla_gemm.
 * A BLAS operation workspace encapsulates the internal state of the operation
 * and provides the execute method, using CBLAS interface.
 * \section tests Tests and examples
 * - \ref examples/cblas_test.c
 * - \ref examples/gemm_block.c
 * \{
 */
#ifdef __cplusplus
extern "C" {
#endif

#ifdef ENABLE_BLAS_SP
  /*! Single-precision floating-point data type */
  typedef float hpnla_float;
  /*! Single-precision floating-point data type for MPI */
#define HDNLA_MPI_FLOAT MPI_FLOAT
#else
  /*! Double-precision floating-point data type */
  typedef double hpnla_float;
  /*! Double-precision floating-point data type for MPI */
#define HDNLA_MPI_FLOAT MPI_DOUBLE
#endif

  /*! GEMM workspace */
  typedef struct hpnla_gemm {
    /*! configuration */
	  hpnla_process_conf* conf;
    /* time recording*/
    float time;
    /*! params */
    void* params;
  } hpnla_gemm;

  /*! Allocate GEMM workspace */
  hpnla_gemm* hpnla_gemm_alloc(hpnla_process_conf* conf);

  /*! Execute GEMM */
  void hpnla_gemm_execute(struct hpnla_gemm* gemm,
      const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
      const int M, const int N, const int K,
      const hpnla_float alpha, const hpnla_float *A, const int lda,
      const hpnla_float *B, const int ldb,  const hpnla_float beta,
      hpnla_float *C, const int ldc);

  /*! Free GEMM workspace */
  void hpnla_gemm_free(hpnla_gemm* gemm);

#ifdef __cplusplus
}
#endif
/*!
 * \}
 */

#endif /* HDNLA_CBLAS_H_ */
