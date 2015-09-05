#include "utils.h"
#include <string.h>

/* variable length args */
#include <stdarg.h>
#include <errno.h>

/* print message to stderr */
void himpi_err(const char *fmt, ...) {
	va_list argp;
	fprintf(stderr, "HiMPI ERROR: rank %d on %s: ", himpi_my_rank_world,
			himpi_my_hostname);
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n");
}

/* print message to stdout if himpi_debug is set and it is >= level */
void himpi_dbg(int level, const char *fmt, ...) {
	va_list argp;
	if (level == 0 || (himpi_debug > 0 && himpi_debug >= level)) {
		fprintf(stdout, "HiMPI: rank %d on %s: ", himpi_my_rank_world,
				himpi_my_hostname);
		va_start(argp, fmt);
		vfprintf(stdout, fmt, argp);
		va_end(argp);
		fprintf(stdout, "\n");
	}
}

/* print abort message and call MPI_Abort to kill run */
void himpi_abort(int rc, const char *fmt, ...) {
	va_list argp;
	fprintf(stderr, "HiMPI ABORT: rank %d on %s: ", himpi_my_rank_world,
			himpi_my_hostname);
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n");

	MPI_Abort(MPI_COMM_WORLD, 0);
}

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

himpi_conf* himpi_get_conf_all(const char* filename, int* num_lines) {
	if (filename == NULL) {
		himpi_abort(-1, "himpi_get_conf_all: filename is null @ %s:%d",
			__FILE__, __LINE__);
	}

	himpi_conf* confs = NULL;
	FILE* stream;
	stream = fopen(filename, "r");
	if (stream == NULL) {
		himpi_err("Can't open the configuration file %s\n", filename);
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

		confs = realloc(confs, sizeof(himpi_conf) * (n + 1));
		int pos = 0;
		int err = sscanf(line, "%d %d %d %d %d %d %d %n", &confs[n].num_procs,
				&confs[n].num_groups, &confs[n].num_levels,
				&confs[n].message_size, &confs[n].alg_in, &confs[n].alg_out,
				&confs[n].op_id, &pos);
		if (err < 6) {
			himpi_abort(-1, "%s: Error reading line: %d \"%s\" err:%d @ %s:%d",
						__func__, __FILE__, __LINE__);
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
 * \return HiMPI default configuration values
 */
himpi_conf default_config(int num_procs, int msg_size) {
	himpi_conf config;
	config.num_procs = num_procs;
	config.num_groups = HIMPI_NUM_GROUPS;
	config.num_levels = HIMPI_NUM_LEVELS;
	config.message_size = msg_size;
	config.alg_in = HIMPI_ALG_IN;
	config.alg_out = HIMPI_ALG_OUT;
	config.op_id = op_all;
	return config;
}

/*!
 * \param confs HiMPI configuration parameters
 * \param num_procs number of processes for which the configuration is being searched
 * \param msg_size the message size for which the configuration is being searched
 * \param op_id HiMPI collective operation id
 * \param num_lines number of lines in the configuration file
 * \param config_found if configuration exists for the given number of process
 * and message size set it to one otherwise set it to zero
 * \return if config found return it otherwise return default config
 */
himpi_conf find_config(himpi_conf* confs, int num_procs, int msg_size,
		himpi_operations op_id, int num_lines, int* config_found) {
	himpi_conf my_conf = default_config(num_procs, msg_size);

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
		int comm_size, himpi_operations op_id, const char* filename) {
	int num_lines;
	himpi_conf* confs = himpi_get_conf_all(filename, &num_lines);
	int p = 0, m = min_msg_size;

	int searched_lines = 0;
	for (m = min_msg_size; m <= max_msg_size; m *= msg_stride) {
		searched_lines++;
	}

	int num_procs = 0;
	for (p = HIMPI_MIN_PROCS; p <= comm_size; p *= 2) {
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
		for (p = comm_size; p >= HIMPI_MIN_PROCS; p /= 2) {
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

himpi_conf himpi_get_my_conf(MPI_Comm comm, int msg_size, int root,
		const char* filename, himpi_operations op_id) {
	int num_lines;
	char* conf_file_name = create_file_name(filename, op_id);
	himpi_conf* confs = himpi_get_conf_all(conf_file_name, &num_lines);

	int num_procs;
	int my_rank;
	MPI_Comm_size(comm, &num_procs);
	MPI_Comm_rank(comm, &my_rank);

	int config_found = 0;
	himpi_conf my_conf = find_config(confs, num_procs, msg_size, op_id,
			num_lines, &config_found);

	if (!config_found) {
		save_himpi_optimal_groups(msg_size, msg_size, henv.msg_stride, root,
				henv.num_levels, henv.bcast_alg_in, henv.bcast_alg_out, op_id,
				1, filename); //TODO: 1
		MPI_Barrier(comm);

		confs = himpi_get_conf_all(conf_file_name, &num_lines);
		my_conf = find_config(confs, num_procs, msg_size, op_id, num_lines,
				&config_found);
	}

	if (!config_found && !my_rank) {
		fprintf(stdout, "No config found. Using bcast. [p: %d, msg: %d]\n",
				num_procs, msg_size);
	}

	return my_conf;
}

void himpi_print_conf(FILE* file, int* num_procs, int* num_groups,
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

