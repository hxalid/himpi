#include "utils.h"

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
		int err = sscanf(line, "%d %d %d %d %d %d %n", &confs[n].num_procs,
				&confs[n].num_groups, &confs[n].num_levels, &confs[n].message_size, &confs[n].alg_in, &confs[n].alg_out, &pos);
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



hmpi_conf hmpi_get_my_conf(MPI_Comm comm, int msg_size, int root, const char* filename, hmpi_operations operation) {
	//get number of groups per process.
	int num_lines;
	hmpi_conf* confs = hmpi_get_conf_all(filename, &num_lines);
	if (confs == NULL) {
		fprintf(stdout, "No config found. Generating config file\n");
		/*
		 * TODO:
		 * at the moment we use only one level of hierarchy,
		 * use algs 0 to ignore selection from our implementations.
		 */
		save_hmpi_optimal_groups(msg_size, root, comm, 1, 0, 0, operation);
		MPI_Barrier(comm);
		confs = hmpi_get_conf_all(filename, &num_lines);
		printf("confs[0].message_size=%d, confs[0].num_groups=%d\n", confs[0].message_size, confs[0].num_groups);
	}


	hmpi_conf my_conf;
	int num_procs;
	MPI_Comm_size(comm, &num_procs);

	int i;
	int config_found = 0;

	my_conf.num_procs=num_procs;
	my_conf.num_groups=1;
	my_conf.num_levels=1;
	my_conf.message_size=msg_size;
	my_conf.alg_in=0;
	my_conf.alg_out=0;

	//TODO
	if (confs != NULL) {
		for(i = 0; i < num_lines; i++){
			if (confs[i].num_procs == num_procs && confs[i].message_size==msg_size) {
				my_conf = confs[i];
				config_found = 1;
				break;
			}
		}
	}

	if (!config_found) {
	   fprintf(stdout, "No config found.\n");
	}

	return my_conf;
}




void hmpi_print_conf(FILE* file, int* num_procs, int* num_groups, int* num_levels, int* msg_sizes,
		int* alg_in, int* alg_out, int size){
	if(file == NULL) file = stdout;

	int i=0;
	fprintf(file, "#num_procs\tnum_groups\tnum_levels\tmsg_size\talg_in\talg_out\n");
	for(i = 0; i < size; i++){
		fprintf(file, "%d\t%d\t%d\t%d\n", num_procs[i], num_groups[i], num_levels[i], msg_sizes[i], alg_in[i], alg_out[i]);
	}

	return;
}






