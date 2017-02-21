#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <mpi.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>

#include "optimal_groups.h"

/*
 * List of algorithms are used in OpenMPI
 */

typedef enum ompi_bcast_algorithm_id {
	DEFAULT, /* 0 */
	BASIC_LINEAR, /* 1 */
	CHAIN, /* 2 */
	PIPELINE, /* 3 */
	SPLIT_BINARY_TREE, /* 4 */
	BINARY_TREE, /* 5 */
	BINOMIAL, /* 6 */
} ompi_bcast_algorithm_id;

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

/* When bcast algorithm id is 0, MPI_Bcast use ompi_coll_tuned_bcast_intra_dec_fixed.
 *  Using this function you can get bcast algorithm which
 *  ompi_coll_tuned_bcast_intra_dec_fixed uses inside.
 *   */
int get_default_bcat_alg_id(int message_size, int communicator_size, int count );

/* Set broadcast algorithm before creation sub communicator */
int set_bcast_alg_t(int alg_id);
