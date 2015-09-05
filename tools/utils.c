#include "utils.h"
#include <string.h>

int validate_groups(int num_groups, int num_procs) {
	if (num_groups <= 1 || num_groups >= num_procs)
		return 0;
	else
		return 1;
}

char *create_rand_elms(int num_elements) {
	char *rand_nums = (char *) malloc(sizeof(char) * num_elements);
	assert(rand_nums != NULL);
	int i;
	for (i = 0; i < num_elements; i++) {
		//TODO
		rand_nums[i] = (rand() / 1024);
	}
	return rand_nums;
}

double hdnla_conf_int(double cl, int reps, double* T) {
	return fabs(gsl_cdf_tdist_Pinv(cl, reps - 1)) * gsl_stats_sd(T, 1, reps)
			/ sqrt(reps);
}

hmpi_conf* hmpi_get_conf_all(const char* filename, int* num_lines) {
	if (filename == NULL) {
		fprintf(stderr, "Error filename null %s\n", __func__);
		MPI_Abort(MPI_COMM_WORLD, 200);
	}

	hmpi_conf* confs = NULL;
	FILE* stream;
	stream = fopen(filename, "r");
	if (stream == NULL) {
		fprintf(stderr, "Can't open the configuration file %s\n", filename);
		//perror("fopen");
		//MPI_Abort(MPI_COMM_WORLD, 201);
		*num_lines = 0;
		return confs;
	}

	char* line = NULL;
	size_t number = 0;
	int err;
	int n = 0;

	while ((err = getline(&line, &number, stream)) != -1) {
		// skip commented lines
		if (line[0] == '#' || line[0] == '\n')
			continue;

		confs = realloc(confs, sizeof(hmpi_conf) * (n + 1));
		int pos = 0;
		int err = sscanf(line, "%d %d %d %d %d %d %d %n", &confs[n].num_procs,
				&confs[n].num_groups, &confs[n].num_levels,
				&confs[n].message_size, &confs[n].alg_in, &confs[n].alg_out,
				&confs[n].op_id, &pos);
		if (err < 6) {
			fprintf(stderr, "Error reading line%d: \"%s\" err:%d in %s\n", n,
					line, err, __func__);
			MPI_Abort(MPI_COMM_WORLD, 204);
		}

		n++;
	}

	fclose(stream);
	free(line);
	*num_lines = n;
	return confs;
}

/*!
 * \param num_procs number of MPI processes
 * \param msg_size the message size in HMPI collective operation
 * \return HMPI default configuration values
 */
hmpi_conf default_config(int num_procs, int msg_size) {
	hmpi_conf config;
	config.num_procs = num_procs;
	config.num_groups = HMPI_NUM_GROUPS;
	config.num_levels = HMPI_NUM_LEVELS;
	config.message_size = msg_size;
	config.alg_in = HMPI_ALG_IN;
	config.alg_out = HMPI_ALG_OUT;
	config.op_id = op_all;
	return config;
}

/*!
 * \param confs hmpi configuration parameters
 * \param num_procs number of processes for which the configuration is being searched
 * \param msg_size the message size for which the configuration is being searched
 * \param op_id HMPI collective operation id
 * \param num_lines number of lines in the configuration file
 * \param config_found if configuration exists for the given number of process
 * and message size set it to one otherwise set it to zero
 * \return if config found return it otherwise return default config
 */
hmpi_conf find_config(hmpi_conf* confs, int num_procs, int msg_size,
		hmpi_operations op_id, int num_lines, int* config_found) {
	hmpi_conf my_conf = default_config(num_procs, msg_size);

	if (confs != NULL) {
		int i;
		for (i = 0; i < num_lines; i++) {
			if ((confs[i].num_procs == num_procs
					&& confs[i].message_size == msg_size)
					&& confs[i].op_id == op_id) {
				my_conf = confs[i];
				*config_found = 1;
				break;
			}
		}
	}
	return my_conf;
}

/*!
 * \param min_msg_size minimum message size in the configuration
 * \param max_msg_size maximum message size in the configuration
 * \param msg_stride the stride to get all message sizes within the min and max
 * \param comm_size MPI communicator size
 * \param op_id HMPI collective operation id
 * \param filename filename to check if requested configuration exists in it
 * \return if the file contains the configuration for the given number of processes
 * and message size then return one otherwise zero
 */
int is_same_config(int min_msg_size, int max_msg_size, int msg_stride,
		int comm_size, hmpi_operations op_id, const char* filename) {
	int num_lines;
	hmpi_conf* confs = hmpi_get_conf_all(filename, &num_lines);
	int p = 0, m = min_msg_size;

	int searched_lines = 0;
	for (m = min_msg_size; m <= max_msg_size; m *= msg_stride) {
		searched_lines++;
	}

	int num_procs = 0;
	for (p = HMPI_MIN_PROCS; p <= comm_size; p *= 2) {
		num_procs++;
	}

	if (num_lines != num_procs * searched_lines) {
		fprintf(stdout,
				"Should generate new config file. num_procs=%d, comm_size: %d, searched_lines: %d num_lines: %d\n",
				num_procs, comm_size, searched_lines, num_lines);
		return 0;
	}

	int matched_lines = 0;
	if (confs != NULL) {
		int i = 0;
		for (p = comm_size; p >= HMPI_MIN_PROCS; p /= 2) {
			for (m = min_msg_size; m <= max_msg_size; m *= msg_stride) {
				if ((confs[i].num_procs == p && confs[i].message_size == m)
						&& confs[i].op_id == op_id) {
					matched_lines++;
				}
				i++;
			}
		}
	}

	if (num_lines != 0 && matched_lines == num_lines) {
		fprintf(stdout, "The same config file [%s] exists...\n", filename);
		return 1;
	} else {
		fprintf(stdout,
				"Should generate new config file. matched_lines: %d, num_lines: %d\n",
				matched_lines, num_lines);
		return 0;
	}
}


hmpi_conf hmpi_get_my_conf(MPI_Comm comm, int msg_size, int root,
		const char* filename, hmpi_operations op_id) {
	int num_lines;
	char* conf_file_name = create_file_name(filename, op_id);
	hmpi_conf* confs = hmpi_get_conf_all(conf_file_name, &num_lines);

	int num_procs;
	int my_rank;
	MPI_Comm_size(comm, &num_procs);
	MPI_Comm_rank(comm, &my_rank);

	int config_found = 0;
	hmpi_conf my_conf = find_config(confs, num_procs, msg_size, op_id,
			num_lines, &config_found);

	if (!config_found) {
		save_hmpi_optimal_groups(msg_size, msg_size, henv.msg_stride, root,
				comm, henv.num_levels, henv.bcast_alg_in, henv.bcast_alg_out,
				op_id, 1, filename); //TODO: 1
		MPI_Barrier(comm);

		confs = hmpi_get_conf_all(conf_file_name, &num_lines);
		my_conf = find_config(confs, num_procs, msg_size, op_id, num_lines,
				&config_found);
	}

	if (!config_found && !my_rank) {
		fprintf(stdout, "No config found. Using bcast. [p: %d, msg: %d]\n",
				num_procs, msg_size);
	}

	return my_conf;
}

void hmpi_print_conf(FILE* file, int* num_procs, int* num_groups,
		int* num_levels, int* msg_sizes, int* alg_in, int* alg_out, int* op_id,
		int size) {
	if (file == NULL)
		file = stdout;

	int i = 0;
	fprintf(file,
			"#num_procs\tnum_groups\tnum_levels\tmsg_size\talg_in\talg_out\top_id\n");
	for (i = 0; i < size; i++) {
		fprintf(file, "%d\t%d\t%d\t%d\n", num_procs[i], num_groups[i],
				num_levels[i], msg_sizes[i], alg_in[i], alg_out[i], op_id[i]);
	}

	return;
}

char* create_file_name(char* file_name, int op_id) {
	const int op_str_length = 3;
	char op_str[op_str_length];
	snprintf(op_str, op_str_length, "-%d", op_id);
	size_t file_length = strlen(file_name) + strlen(op_str) + 1;
	char* conf_file_name = (char*) malloc(file_length);
	assert(conf_file_name != NULL);
	snprintf(conf_file_name, file_length, "%s%s", file_name, op_str);
	return conf_file_name;
}

