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
	FILE* stream;
	stream = fopen(filename, "r");
	if (stream == NULL) {
		printf("Try to open the configuration file %s\n", filename);
		perror("fopen");
		MPI_Abort(MPI_COMM_WORLD, 201);
	}

	char* line = NULL;
	size_t number = 0;
	int err;
	int n = 0;
	hmpi_conf* confs = NULL;

	while ((err = getline(&line, &number, stream)) != -1) {
		// skip commented lines
		if (line[0] == '#' || line[0] == '\n')
			continue;

		confs = realloc(confs, sizeof(hmpi_conf) * (n + 1));
		int pos = 0;
		int err = sscanf(line, "%d %d %d %d %d %n", &confs[n].num_procs,
				&confs[n].num_groups, &confs[n].num_levels, &confs[n].alg_in, &confs[n].alg_out, &pos);
		if (err < 5) {
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



hmpi_conf hmpi_get_my_conf(MPI_Comm comm, const char* filename) {
	//get number of groups per process.
	int num_lines;
	hmpi_conf* confs = hmpi_get_conf_all(filename, &num_lines);
	hmpi_conf my_conf;
	int num_procs;
	MPI_Comm_size(comm, &num_procs);

	int i;

	//TODO
	if (confs == NULL) {
		my_conf.num_procs=-1;
		my_conf.num_groups=-1;
		my_conf.num_levels=-1;
		my_conf.alg_in=0;
		my_conf.alg_out=0;
	}


	for(i = 0; i < num_lines; i++){
		if (confs[i].num_procs == num_procs) {
			my_conf = confs[i];
			break;
		} else {
			my_conf.num_procs=-1;
			my_conf.num_groups=-1;
			my_conf.num_levels=-1;
			my_conf.alg_in=0;
			my_conf.alg_out=0;
		}
	}

	return my_conf;
}




void hmpi_print_conf(FILE* file, int* num_procs, int* num_groups, int* num_levels,
		int* alg_in, int* alg_out, int size){
	if(file == NULL) file = stdout;

	int i=0;
	fprintf(file, "#num_procs\tnum_groups\tnum_levels\talg_in\talg_out\n");
	for(i = 0; i < size; i++){
		fprintf(file, "%d\t%d\t%d\t%d\n", num_procs[i], num_groups[i], num_levels[i], alg_in[i], alg_out[i]);
	}

	return;
}






