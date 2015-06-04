/*
 * hmpi_init.c
 *
 *  Created on: 29 Apr 2015
 *      Author: khalidhasanov
 */

#include <stdlib.h>
#include "hmpi.h"
#include "tools/optimal_groups.h"

const char *HMPI_CONF_FILE_NAME;
hmpi_env henv;

void init_hmpi_env(void) {
	henv.bcast_alg_in =
			(getenv("HMPI_BCAST_ALG_IN") != NULL) ?
					strtol(getenv("HMPI_BCAST_ALG_IN"), NULL, 10) :
					HMPI_ALG_IN;
	henv.bcast_alg_out =
			(getenv("HMPI_BCAST_ALG_OUT") != NULL) ?
					strtol(getenv("HMPI_BCAST_ALG_OUT"), NULL, 10) :
					HMPI_ALG_OUT;
	henv.min_msg =
			(getenv("HMPI_MIN_MSG") != NULL) ?
					strtol(getenv("HMPI_MIN_MSG"), NULL, 10) : HMPI_MIN_MSG;
	henv.max_msg =
			(getenv("HMPI_MAX_MSG") != NULL) ?
					strtol(getenv("HMPI_MAX_MSG"), NULL, 10) : HMPI_MAX_MSG;
	henv.msg_stride =
			(getenv("HMPI_MSG_STRIDE") != NULL) ?
					strtol(getenv("HMPI_MSG_STRIDE"), NULL, 10) :
					HMPI_MSG_STRIDE;
	henv.root =
			(getenv("HMPI_ROOT_PROC") != NULL) ?
					strtol(getenv("HMPI_ROOT_PROC"), NULL, 10) : HMPI_ROOT_PROC;
	henv.num_levels =
			(getenv("HMPI_NUM_LEVELS") != NULL) ?
					strtol(getenv("HMPI_NUM_LEVELS"), NULL, 10) :
					HMPI_NUM_LEVELS;
	henv.generate_config =
			(getenv("HMPI_GENERATE_CONFIG") != NULL) ?
					strtol(getenv("HMPI_GENERATE_CONFIG"), NULL, 10) :
					HMPI_GENERATE_CONFIG;

	henv.conf_file_name = getenv("HMPI_CONF_FILE");
}

int HMPI_Init(int *argc, char ***argv) {
	init_hmpi_env();

	int rc = MPI_Init(argc, argv);

	if (rc != MPI_SUCCESS) {
		return rc;
	}

#ifdef HMPI_GROUP_CONFIG
	HMPI_CONF_FILE_NAME = HMPI_GROUP_CONFIG;
#else
	HMPI_CONF_FILE_NAME = HMPI_CONF_FILE;
	//fprintf(stdout, "No group config file specified, g=sqrt(p) will be used\n");
#endif

	/*
	 * User can overwrite the config file name by environment variable
	 */
	if (henv.conf_file_name!=NULL)
		HMPI_CONF_FILE_NAME = henv.conf_file_name;

	int comm_size;
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	int should_generate_config = 0;
#ifdef HMPI_GENERATE_CONFIG
	should_generate_config = 1;
#endif

	if (henv.generate_config == 0 || henv.generate_config == 1)
		should_generate_config = henv.generate_config;

	if (should_generate_config)
		save_hmpi_optimal_groups(henv.min_msg, henv.max_msg, henv.msg_stride,
				henv.root, MPI_COMM_WORLD, henv.num_levels, henv.bcast_alg_in,
				henv.bcast_alg_out, op_bcast, 0); //TODO: 0

	return rc;
}
