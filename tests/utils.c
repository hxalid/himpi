#include "utils.h"

int validate_input(int num_groups, int num_procs) {
    if (num_groups <= 1 || num_groups >= num_procs)
        return 0;
    else
        return 1;
}


char *create_rand_elms(int num_elements) {
    char *rand_nums = (char *) malloc(sizeof (char) * num_elements);
    assert(rand_nums != NULL);
    int i;
    for (i = 0; i < num_elements; i++) {
        rand_nums[i] = (rand() / (int) MAX_INT_SIZE);
    }
    return rand_nums;
}

