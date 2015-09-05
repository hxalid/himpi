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
					strtol(getenv("HIMPI_ROOT_PROC"), NULL, 10) : HIMPI_ROOT_PROC;
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

#ifdef HiMPI_GROUP_CONFIG
	himpi_conf_file_name = HiMPI_GROUP_CONFIG;
#else
	himpi_conf_file_name = HIMPI_CONF_FILE;
	//fprintf(stdout, "No group config file specified, g=sqrt(p) will be used\n");
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

	int is_env_value = 0;
	if (getenv("HIMPI_OPID") != NULL) {
		himpi_operation = strtol(getenv("HIMPI_OPID"), NULL, 10);
	}

	if (himpi_operation < op_bcast || himpi_operation > op_all) {
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank == 0)
			fprintf(stdout,
					"Wrong HiMPI operation id: %d given via configuration or HiMPI_OPID environment. "
					"Using %d instead\n"
					"Allowed operation ids: [0: bcast, 1: reduce, 2: allreduce, "
					"3: scatter, 4: gather, 5: all the previous collectives]\n", himpi_operation, op_all);
		himpi_operation = op_all;
		//MPI_Abort(MPI_COMM_WORLD, -1); //TODO: -1
	}

	int comm_size;
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//TODO: this part should be enabled/disabled by some debug level
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);

	fprintf(stdout, "[%d, %s] -> fayil=%s\n", rank, processor_name,
			himpi_conf_file_name);

	//

	int should_generate_config = 0;
#ifdef HMPI_GENERATE_CONFIG
	should_generate_config = 1;
#endif

	if (getenv("HIMPI_GENERATE_CONFIG") != NULL) {
		should_generate_config = henv.generate_config;
	}


	if (should_generate_config) {
		if (himpi_operation != op_all) {
			save_himpi_optimal_groups(henv.min_msg, henv.max_msg,
					henv.msg_stride, henv.root, MPI_COMM_WORLD, henv.num_levels,
					henv.bcast_alg_in, henv.bcast_alg_out, himpi_operation, 0, himpi_conf_file_name); //TODO: 0
		}
		else {
			int op = op_bcast;
			for (; op < op_all; op++) {
				save_himpi_optimal_groups(henv.min_msg, henv.max_msg,
						henv.msg_stride, henv.root, MPI_COMM_WORLD,
						henv.num_levels, henv.bcast_alg_in, henv.bcast_alg_out,
						op, 0, himpi_conf_file_name); //TODO: 0
			}

		}

	}

	return rc;
}
