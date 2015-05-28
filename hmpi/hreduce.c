/*
 * hreduce.c
 *
 *  Created on: 25 Mar 2015
 *      Author: Khalid Hasanov
 */

#include <mpi.h>
#include <stdlib.h>

#include "hmpi.h"
#include "tools/utils.h"

/*
 * TODO:  It works only for root=0. Fix it to work with any possible root.
 */
int hierarchical_reduce(void *snd_buffer, void* rcv_buffer, int count, MPI_Datatype datatype, MPI_Op op,
        int root, MPI_Comm comm_world, int num_groups, int num_levels, int alg_in, int alg_out) {
	int err;
    int pg;
    int rank;
    int size;
    int stride;
    int root_inside;
    int root_outside;
    int my_group;
    double *reduce_out;
    MPI_Comm in_group_comm, out_group_comm;

    MPI_Comm_rank(comm_world, &rank);
    MPI_Comm_size(comm_world, &size);

    if (size == 1) return MPI_SUCCESS;

    if (0 == count) {
        return MPI_SUCCESS;
    }


    /*TODO make num_groups configurable*/
    if (size > HREDUCE_MIN_PROCS && ( num_groups > 1 && num_groups < size ) ) {
        pg = size / num_groups;
        my_group = rank / pg;
        //stride = root / pg;

        MPI_Comm_split(comm_world, (rank - my_group * pg == 0) ? 0 : MPI_UNDEFINED, rank, &out_group_comm);
        MPI_Comm_split(comm_world, my_group, rank, &in_group_comm);


        root_inside = root; // root_inside = stride;
        /*
         * Start reduce inside groups
         */
        switch (num_levels) {
            case 1: // 1 level of hierarchy
                reduce_out = (double*) malloc(count*sizeof(double));

                err = MPI_Reduce(snd_buffer, reduce_out, count, datatype, op, root_inside, in_group_comm);
                break;
            case -1:
                // Just to see if reduce inside groups is better than that of with MPI_COMM_WORLD
                err = MPI_Reduce(snd_buffer, rcv_buffer, count, datatype, op, root_inside, in_group_comm);
                break;
            default: // e.g. -2
                // Don't reduce inside groups. Just to see if reduce between groups is better than that of with MPI_COMM_WORLD
                break;
        }


         /*
         * Start reduce between groups, if num_levels == -1 then reduce only inside
         */
        if (out_group_comm != MPI_COMM_NULL && err == MPI_SUCCESS && num_levels != -1) {
            root_outside = root; //  root / pg;
            err = MPI_Reduce( (num_levels==1)?reduce_out:snd_buffer,  rcv_buffer, count, datatype, op, root_outside, out_group_comm);
        //  MPI_Reduce(&reduce_in, &reduce_out, 1, datatype, op, root_outside, out_group_comm);
        }

        //TODO: check
        if (num_levels==1) free(reduce_out);

        if (out_group_comm != MPI_COMM_NULL)
            MPI_Comm_free(&out_group_comm);
        if (in_group_comm != MPI_COMM_NULL)
            MPI_Comm_free(&in_group_comm);
    } else if (size <= HREDUCE_MIN_PROCS || ( num_groups == 1 || num_groups == size) ) {
        //fprintf(stdout, "Using non-hierarchical reduce\n");
        err = MPI_Reduce(snd_buffer, rcv_buffer, count, datatype, op, root, comm_world);
    } else {
        /*TODO*/
        //fprintf(stdout, "Wrong number of groups\n");
    }

    return err;
}


//TODO
int HMPI_Reduce(void *snd_buffer, void* rcv_buffer, int count, MPI_Datatype datatype, MPI_Op op,
        int root, MPI_Comm comm) {
	MPI_Aint extent, lb;
	MPI_Type_get_extent(datatype, &lb, &extent);
	int msg_size = extent * count;

	/*
	 * TODO: Are you sure all processes have to open config file for read?
	 */
	hmpi_conf my_conf = hmpi_get_my_conf(comm, msg_size, root, HMPI_CONF_FILE_NAME, op_reduce);

	return hierarchical_reduce(snd_buffer, rcv_buffer, count, datatype, op, root, comm, my_conf.num_groups,
			my_conf.num_levels, my_conf.alg_in, my_conf.alg_out);
}









