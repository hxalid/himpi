#ifndef MPIB_SHARED_H_
#define MPIB_SHARED_H_

#include <mpi.h>

/*!
 * \defgroup shared API for shared libraries
 * \{
 */
#ifdef __cplusplus
extern "C" {
#endif

/*! Name of the initialization function in the shared library */
#define MPIB_SHARED_INITIALIZE "mpib_shared_initialize"

/*!
 * Type of functions that parse suboptions and allocate global variables at all processors.
 * Implement the function \c mpib_shared_initialize of this type in the shared library (optional).
 * \param comm MPI communicator
 * \param subopts suboptions
 */
typedef int (*MPIB_shared_initialize)(MPI_Comm comm, char* subopts);

/*! Name of the finalization function in the shared library */
#define MPIB_SHARED_FINALIZE "mpib_shared_finalize"

/*!
 * Type of functions that free global variables at all processors.
 * Implement the function \c mpib_shared_finalize of this type in the shared library (optional).
 * \param comm MPI communicator
 */
typedef int (*MPIB_shared_finalize)(MPI_Comm comm);

#ifdef __cplusplus
}
#endif
/*!
 * \}
 */

#endif /* MPIB_SHARED_H_ */
