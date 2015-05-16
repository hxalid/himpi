#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <mpi.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>


/*! hmpi configuration parameters */
typedef struct hmpi_conf {
	int num_procs;
	int num_groups;
	int num_levels;
	int alg_in;
	int alg_out;
} hmpi_conf;



int validate_groups(int num_groups, int num_procs);

char *create_rand_elms(int num_elements);

double hdnla_conf_int(double cl, int reps, double* T);

hmpi_conf hmpi_get_my_conf(MPI_Comm comm, const char* filename);

void hmpi_print_conf(FILE* file, int* num_procs, int* num_groups, int* num_levels, int* alg_in, int* alg_out, int size);

hmpi_conf* hmpi_get_conf_all(const char* filename, int* num_lines);


