/*!
 * hmpi_init.c
 *
 *  Created on: 29 Apr 2015
 *      Author: Khalid Hasanov
 */

#include <stdlib.h>
#include "himpi.h"
#include "tools/optimal_groups.h"
#include "tools/utils.h"

const char *himpi_conf_file_name;
himpi_env henv;
int himpi_debug = 0;
int himpi_my_rank_world = MPI_PROC_NULL;
int himpi_num_ranks_world;
char himpi_my_hostname[256] = "";
MPI_Comm himpi_comm_world = MPI_COMM_NULL;

void init_himpi_env(void) {
	henv.bcast_alg_in =
			(getenv("HIMPI_BCAST_ALG_IN") != NULL) ?
					strtol(getenv("HIMPI_BCAST_ALG_IN"), NULL, 10) :
					HIMPI_ALG_IN;
	henv.bcast_alg_out =
			(getenv("HIMPI_BCAST_ALG_OUT") != NULL) ?
					strtol(getenv("HIMPI_BCAST_ALG_OUT"), NULL, 10) :
					HIMPI_ALG_OUT;
	henv.min_msg =
			(getenv("HIMPI_MIN_MSG") != NULL) ?
					strtol(getenv("HIMPI_MIN_MSG"), NULL, 10) : HIMPI_MIN_MSG;
	henv.max_msg =
			(getenv("HIMPI_MAX_MSG") != NULL) ?
					strtol(getenv("HIMPI_MAX_MSG"), NULL, 10) : HIMPI_MAX_MSG;
	henv.msg_stride =
			(getenv("HIMPI_MSG_STRIDE") != NULL) ?
					strtol(getenv("HIMPI_MSG_STRIDE"), NULL, 10) :
					HIMPI_MSG_STRIDE;
	henv.root =
			(getenv("HIMPI_ROOT_PROC") != NULL) ?
					strtol(getenv("HIMPI_ROOT_PROC"), NULL, 10) :
					HIMPI_ROOT_PROC;
	henv.num_levels =
			(getenv("HIMPI_NUM_LEVELS") != NULL) ?
					strtol(getenv("HIMPI_NUM_LEVELS"), NULL, 10) :
					HIMPI_NUM_LEVELS;
	henv.generate_config =
			(getenv("HIMPI_GENERATE_CONFIG") != NULL) ?
					strtol(getenv("HIMPI_GENERATE_CONFIG"), NULL, 10) :
					HIMPI_GENERATE_CONFIG;

	henv.conf_file_name = getenv("HIMPI_CONF_FILE");
}

int HiMPI_Init(int *argc, char ***argv) {
	init_himpi_env();

	int rc = MPI_Init(argc, argv);

	if (rc != MPI_SUCCESS) {
		return rc;
	}

	MPI_Comm_dup(MPI_COMM_WORLD, &himpi_comm_world);
	if (himpi_comm_world == MPI_COMM_NULL) {
		himpi_abort(-1,
				"himpi_init: Failed to duplicate MPI_COMM_WORLD @ %s:%d",
				__FILE__, __LINE__);
	}

	MPI_Comm_size(himpi_comm_world, &himpi_num_ranks_world);
	MPI_Comm_rank(himpi_comm_world, &himpi_my_rank_world);

	/* get my hostname */
	if (gethostname(himpi_my_hostname, sizeof(himpi_my_hostname)) != 0) {
		himpi_abort(-1, "himpi_init: Failed to get hostname @ %s:%d",
		__FILE__, __LINE__);
	}

#ifdef HIMPI_GROUP_CONFIG
	himpi_conf_file_name = HIMPI_GROUP_CONFIG;
#else
	himpi_conf_file_name = HIMPI_CONF_FILE;
#endif

	/*!
	 * User can overwrite the config file name by environment variable
	 */
	if (henv.conf_file_name != NULL)
		himpi_conf_file_name = henv.conf_file_name;

	int himpi_operation = op_all;
#ifdef HIMPI_OPID
	himpi_operation = strtol(HIMPI_OPID, NULL, 10);
#endif

	if (getenv("HIMPI_DEBUG") != NULL) {
		himpi_debug = strtol(getenv("HIMPI_DEBUG"), NULL, 10);
	}

	if (getenv("HIMPI_OPID") != NULL) {
		himpi_operation = strtol(getenv("HIMPI_OPID"), NULL, 10);
	}

	if (himpi_operation < op_bcast || himpi_operation > op_all) {
		if (himpi_my_rank_world == 0)
			himpi_dbg(0,
					"Wrong HiMPI operation id: %d given via configuration or HiMPI_OPID environment. "
							"Using %d instead\n"
							"Allowed operation ids: [0: bcast, 1: reduce, 2: allreduce, "
							"3: scatter, 4: gather, 5: all the previous collectives]\n",
					himpi_operation, op_all);
		himpi_operation = op_all;
		//MPI_Abort(MPI_COMM_WORLD, -1); //TODO: -1
	}

	int generate_config = 0;
#ifdef HMPI_GENERATE_CONFIG
	generate_config = 1;
#endif

	if (getenv("HIMPI_GENERATE_CONFIG") != NULL) {
		generate_config = henv.generate_config;
	}

	himpi_dbg(1,
			"HiMPI_Init: initialization is continuing with a call to save_himpi_optimal_groups. "
					"himpi_operation: %d, rank: %d (%d), hostname: %s",
			himpi_operation, himpi_my_rank_world, himpi_num_ranks_world,
			himpi_my_hostname);

	if (generate_config) {
		if (himpi_operation != op_all) {
			save_himpi_optimal_groups(henv.min_msg, henv.max_msg,
					henv.msg_stride, henv.root, henv.num_levels,
					henv.bcast_alg_in, henv.bcast_alg_out, himpi_operation, 0,
					himpi_conf_file_name); //TODO: 0
		} else {
			int op = op_bcast;
			for (; op < op_all; op++) {
				save_himpi_optimal_groups(henv.min_msg, henv.max_msg,
						henv.msg_stride, henv.root, henv.num_levels,
						henv.bcast_alg_in, henv.bcast_alg_out, op, 0,
						himpi_conf_file_name); //TODO: 0
			}

		}

	}

	return rc;
}
