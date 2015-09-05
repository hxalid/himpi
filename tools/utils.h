#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <mpi.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>

#include "optimal_groups.h"

int validate_groups(int num_groups, int num_procs);

char *create_rand_elms(int num_elements);

double hdnla_conf_int(double cl, int reps, double* T);

void himpi_print_conf(FILE* file, int* num_procs, int* num_groups,
		int* num_levels, int* msg_sizes, int* alg_in, int* alg_out, int* op_id, int size);

himpi_conf himpi_get_my_conf(MPI_Comm comm, int msg_size, int root,
		const char* filename, himpi_operations op_id);

himpi_conf* himpi_get_conf_all(const char* filename, int* num_lines);

char* create_file_name(char* file_name, int op_id);

int is_same_config(int min_msg_size, int max_msg_size, int msg_stride,
		int comm_size, himpi_operations op_id, const char* filename);

void himpi_err(const char *fmt, ...);
void himpi_dbg(int level, const char *fmt, ...);
void himpi_abort(int rc, const char *fmt, ...);
