#ifndef MPIB_P2P_BENCHMARK_H_
#define MPIB_P2P_BENCHMARK_H_

#include "mpib_measurement.h"
#include "mpib_defs.h"

/*!
 * \defgroup p2p_benchmarks Point-to-point benchmarks
 * This module provides the point-to-point benchmarks.
 * \{
 */
#ifdef __cplusplus
extern "C" {
#endif

/*!
 * Point-to-point benchmark. Estimates the execution time of the point-to-point communications
 * between a pair of processors in the MPI communicator.
 * Performs series of communication experiments to obtain reliable results.
 * \param send send
 * \param recv recv
 * \param comm communicator, number of nodes should be \f$ \ge 2 \f$
 * \param measure measure processor
 * \param mirror mirror processor
 * \param M message size
 * \param precision measurement precision
 * \param result measurement result (significant only at the measure processor)
 */
void MPIB_measure_p2p(MPIB_Send send, MPIB_Recv recv, MPI_Comm comm, int measure, int mirror,
		int M, MPIB_precision precision, MPIB_result* result);

/*!
 * Point-to-point benchmark. Estimates the execution time of the point-to-point communication
 * between a pair of processors in the MPI communicator for different message sizes.
 * Performs series of communication experiments to obtain reliable results.
 * \param container communication operation container
 * \param comm communicator, number of nodes should be \f$ \ge 2 \f$
 * \param measure measure processor
 * \param mirror mirror processor
 * \param msgset message sizes
 * \param precision measurement precision
 * \param count the number of measurements performed (significant only at the measure processor)
 * \param results array of measurement results (significant only at the measure processor,
 * allocated by this function and must be deallocated by user)
 */
void MPIB_measure_p2p_msgset(MPIB_Send send, MPIB_Recv recv, MPI_Comm comm, int measure, int mirror,
		MPIB_msgset msgset, MPIB_precision precision, int* count, MPIB_result** results);

/*!
 * Parallel point-to-point benchmark for two partitions of processes. Estimates the execution time of the parallel point-to-point communication
 * between two partitions of processors in the MPI communicator for different message sizes.
 * Performs series of communication experiments to obtain reliable results.
 * \param container communication operation container
 * \param comm communicator, number of nodes should be \f$ \ge 2 \f$
 * \param msgset message sizes
 * \param precision measurement precision
 * \param count the number of measurements performed (significant only at the measure processor)
 * \param results array of measurement results (significant only at the measure processor,
 * allocated by this function and must be deallocated by user)
 */
void MPIB_measure_p2p_parallel_msgset(MPIB_Send send, MPIB_Recv recv, MPI_Comm comm,
		MPIB_msgset msgset, MPIB_precision precision, int* count, MPIB_result** results);

/*! \f$ C_n^2 \f$ */
#define MPIB_C2(n) (n) * ((n) - 1) / 2

/*!
 * For a symmetric square matrix stored in the array of \f$ C_n^2 \f$ elements,
 * returns the index of the \f$ (i, j) \f$ element, \f$ i \ne j < n \f$:
 * \f$ \displaystyle\frac{(n - 1) + (n - I)}{2} I + (J - I - 1) \f$, \f$ I = min(i, j) \f$, \f$ J = max(i, j) \f$
 */
#define MPIB_IJ2INDEX(n, i, j) (2 * (n) - ((i) < (j) ? (i) : (j)) - 1) * (((i) < (j) ? (i) : (j))) / 2 + (((i) < (j) ? (j) : (i))) - (((i) < (j) ? (i) : (j))) - 1

/*!
 * Point-to-point benchmark. Estimates the execution time of the point-to-point communications
 * between all pairs of processors in the MPI communicator.
 * Performs series of communication experiments to obtain reliable results.
 * \param container communication operation container
 * \param comm communicator, number of nodes should be \f$ \ge 2 \f$
 * \param parallel several non-overlapped point-to-point communications at the same time if non-zero
 * \param M message size
 * \param precision measurement precision
 * \param results array of \f$ C_n^2 \f$ measurement results
 */
void MPIB_measure_allp2p(MPIB_Send send, MPIB_Recv recv, MPI_Comm comm, int parallel,
		int M, MPIB_precision precision, MPIB_result* results);

#ifdef __cplusplus
}
#endif
/*!
 * \}
 */

#endif /*MPIB_P2P_BENCHMARK_H_*/
