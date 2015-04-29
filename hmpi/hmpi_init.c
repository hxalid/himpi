/*
 * hmpi_init.c
 *
 *  Created on: 29 Apr 2015
 *      Author: khalidhasanov
 */

#include "hmpi.h"


int HMPI_Init( int *argc, char ***argv ) {
	int rc = MPI_Init(argc, argv);

	if (rc != MPI_SUCCESS) {
		return rc;
	}

    //fprintf(stdout, "HMPI_Init finished successfully\n");

    return rc;
}
