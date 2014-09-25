#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>

#define MAX_INT_SIZE 100

int validate_input(int num_groups, int num_procs);

char *create_rand_elms(int num_elements);

double hdnla_conf_int(double cl, int reps, double* T);
