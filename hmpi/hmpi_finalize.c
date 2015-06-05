/*
 * hmpi_finalize.c
 *
 *  Created on: 30 Apr 2015
 *      Author: Khalid Hasanov
 */

#include "hmpi.h"

int HMPI_Finalize( void ) {

	return MPI_Finalize();
}
