#ifndef MPIB_COLL_CONTAINERS_H_
#define MPIB_COLL_CONTAINERS_H_

#include "mpib_coll_benchmarks.h"
#include "mpib_defs.h"

/*!
 * \defgroup coll_containers Containers for collective communication operations
 * \{
 */
#ifdef __cplusplus
extern "C" {
#endif

/*! Frees collective container */
void MPIB_coll_container_free(MPIB_coll_container* container);

/*! Allocates Scatter container */
MPIB_coll_container* MPIB_Scatter_container_alloc(MPIB_Scatter scatter);

/*! Allocates Gather container */
MPIB_coll_container* MPIB_Gather_container_alloc(MPIB_Gather gather);

/*! Allocates Bcast container */
MPIB_coll_container* MPIB_Bcast_container_alloc(MPIB_Bcast bcast);

/*! Allocates Reduce container */
MPIB_coll_container* MPIB_Reduce_container_alloc(MPIB_Reduce reduce);

/*! Allocates MPI_Comm_dup-MPI_Comm_free container */
MPIB_coll_container* MPIB_Comm_dup_free_container_alloc();

/*! Allocates Scatterv container */
MPIB_coll_container* MPIB_Scatterv_container_alloc(MPIB_Scatterv scatterv,
		const double* factors);

/*! Allocates Gatherv container */
MPIB_coll_container* MPIB_Gatherv_container_alloc(MPIB_Gatherv gatherv,
		const double* factors);

/*! Allocates Gatherv container */
MPIB_coll_container* MPIB_Alltoall_container_alloc(MPIB_Alltoall alltoall);

/*HMPI_Bcast container*/
MPIB_coll_container* MPIB_HBcast_container_alloc(MPIB_HBcast hbcast,
		int num_groups, int num_levels, int alg_in, int alg_out);

MPIB_coll_container* MPIB_HReduce_container_alloc(MPIB_HReduce hreduce,
		int num_groups, int num_levels, int alg_in, int alg_out);

MPIB_coll_container* MPIB_HAllreduce_container_alloc(MPIB_HAllreduce hallreduce,
		int num_groups, int num_levels, int alg_in, int alg_out);

MPIB_coll_container* MPIB_HGather_container_alloc(MPIB_HGather hgather,
		int num_groups, int num_levels, int alg_in, int alg_out);

MPIB_coll_container* MPIB_HScatter_container_alloc(MPIB_HScatter hscatter,
		int num_groups, int num_levels, int alg_in, int alg_out);

#ifdef __cplusplus
}
#endif
/*!
 * \}
 */

#endif /* MPIB_COLL_CONTAINERS_H_ */
