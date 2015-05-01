/*
 * hmpi_init.c
 *
 *  Created on: 29 Apr 2015
 *      Author: khalidhasanov
 */

#include "hmpi.h"

const char *HMPI_CONF_FILE_NAME = "fayil.conf";

/*
 * num_groups is an extern variable defined in hmpi.h
 */
int num_groups;

int HMPI_Init( int *argc, char ***argv ) {
	int rc = MPI_Init(argc, argv);

	if (rc != MPI_SUCCESS) {
		return rc;
	}


#ifdef HMPI_GROUP_CONFIG
	printf("HMPI_GROUP_CONFIG=%s\n", HMPI_GROUP_CONFIG);
	num_groups = hmpi_get_num_groups(MPI_COMM_WORLD, HMPI_CONF_FILE_NAME);
#else
	num_groups = 0;
#endif

	printf("num_groups=%d\n", num_groups);

	//Find out optimal number of groups here

	//int comm_size;
	//MPI_Comm_size(MPI_COMM_WORLD, &comm_size);


	//HMPI_Bcast(array, msg_size, MPI_CHAR, root, MPI_COMM_WORLD, rec, alg);

    //fprintf(stdout, "HMPI_Init finished successfully. comm_size was: %d\n", comm_size);

    return rc;
}
