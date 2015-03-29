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
		rand_nums[i] = (rand() / (int) MAX_INT_SIZE);
	}
	return rand_nums;
}

double hdnla_conf_int(double cl, int reps, double* T) {
	return fabs(gsl_cdf_tdist_Pinv(cl, reps - 1)) * gsl_stats_sd(T, 1, reps)
			/ sqrt(reps);
}

hmpi_conf* hmpi_get_conf_all(char* filename, int* size) {
	if (filename == NULL) {
		fprintf(stderr, "Error filename null %s\n", __func__);
		MPI_Abort(MPI_COMM_WORLD, 200);
	}
	FILE* stream;
	stream = fopen(filename, "r");
	if (stream == NULL) {
		debug_print(0, (size_t) 0, (size_t) 0,
				"Try to open the configuration file %s\n", filename);
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
		int err = sscanf(line, "%s %d %s %s %s %n", confs[n].num_procs,
				&confs[n].num_groups, confs[n].num_levels, &pos);
		if (err < 4) {
			fprintf(stderr, "Error reading line%d: \"%s\" err:%d in %s\n", n,
					line, err, __func__);
			MPI_Abort(MPI_COMM_WORLD, 204);
		}

		n++;
	}

	free(line);
	*size = n;
	return confs;
}



void hmpi_print_conf(FILE* file, int* num_procs, int* num_groups, int* num_levels, int size){
	if(file == NULL) file = stdout;

	int i=0;
	fprintf(file, "#num_procs\tnum_groups\tnum_levels\n");
	for(i = 0; i < size; i++){
		fprintf(file, "%d\t%d\t%d\n", num_procs[i], num_groups[i], num_levels[i]);
	}

	return;
}






