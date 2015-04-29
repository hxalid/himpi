/*
 * hmpi_finalize.c
 *
 *  Created on: 30 Apr 2015
 *      Author: khalidhasanov
 */

#include "hmpi.h"

int HMPI_Finalize( void ) {

	//fprintf(stdout, "Finished\n");

	return MPI_Finalize();
}
