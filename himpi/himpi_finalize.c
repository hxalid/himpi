/*
 * hmpi_finalize.c
 *
 *  Created on: 30 Apr 2015
 *      Author: Khalid Hasanov
 */

#include "himpi.h"

int HiMPI_Finalize( void ) {
	if (himpi_comm_world != MPI_COMM_NULL)
		MPI_Comm_free(&himpi_comm_world);
	return MPI_Finalize();
}
