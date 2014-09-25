#include "hbcast.h"
#include "utils.h"


/*
 * TODO:  It works only for root=0. Fix it to work with any possible root.
 */
t_bcast_response MPI_HBcast(void *buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm_world, int num_groups, int rec, int alg, int debug) {
    int pg;
    int rank;
    int size;
    int stride;
    int root_inside;
    int root_outside;
    int my_group;
    t_bcast_response bcast_response = {
        .rec_out_group_comm = -1,
        .rec_in_group_comm = -1,
        .rec_comm_world = -1,
        .err_flag = 0
    };

    MPI_Comm in_group_comm, out_group_comm;

    MPI_Comm_rank(comm_world, &rank);
    MPI_Comm_size(comm_world, &size);
    if (size == 1) return bcast_response;

    /*TODO make num_groups configurable*/
    if (size > HBCAST_MIN_PROCS && validate_input(num_groups, size)) {
        pg = size / num_groups;
        my_group = rank / pg;
        stride = root / pg;

        MPI_Comm_split(comm_world, (rank - my_group * pg == 0) ? 0 : MPI_UNDEFINED, rank, &out_group_comm);
        //  MPI_Comm_split(comm_world, (rank - my_group * pg == stride) ? 0 : MPI_UNDEFINED, rank, &out_group_comm);
        MPI_Comm_split(comm_world, my_group, rank, &in_group_comm);

       /* if (debug == 2)
           MPIX_Get_property(in_group_comm, MPIDO_RECT_COMM, &(bcast_response.rec_in_group_comm));
        */

        /*
         * Start broadcast between groups
         */
        if (out_group_comm != MPI_COMM_NULL && rec != -1) { // if rec == -1 then broadcast only inside
            int out_size=0;
            int out_rank=-1;
            root_outside = root; //  root / pg; 
            hpnla_bcast(buffer, count, datatype, root_outside, out_group_comm, alg);
            if (num_groups==-99) {
            MPI_Comm_size(out_group_comm, &out_size);
            MPI_Comm_rank(out_group_comm, &out_rank);
            printf("out_rank=%d, out_size=%d\n", out_rank, out_size);
            }
           /* if (debug == 2)
               MPIX_Get_property(out_group_comm, MPIDO_RECT_COMM, &(bcast_response.rec_out_group_comm));
            */ 
        }

        root_inside = root; // root_inside = stride;
        /*
         * Start broadcast inside groups
         */
        switch (rec) {
            case 3: // more levels of hierarchy, until size < HBCAST_MIN_PROCS
                MPI_HBcast(buffer, count, datatype, root_inside, in_group_comm, num_groups, 3, alg, debug);
                break;
            case 2: // two levels of hierarchy. Try making hierarchy further inside groups
                MPI_HBcast(buffer, count, datatype, root_inside, in_group_comm, num_groups, 1, alg, debug);
                break;
            case 1: // 1 level of hierarchy
                hpnla_bcast(buffer, count, datatype, root_inside, in_group_comm, alg);
                break;
            case -1:
                // Just to see if broadcast inside groups is better than that of with MPI_COMM_WORLD
                hpnla_bcast(buffer, count, datatype, root_inside, in_group_comm, alg);
                break;
            default: // e.g. -2
                // Don't broadcast inside groups. Just to see if broadcast between groups is better than that of with MPI_COMM_WORLD
                break;
        }

        if (out_group_comm != MPI_COMM_NULL)
            MPI_Comm_free(&out_group_comm);
        if (in_group_comm != MPI_COMM_NULL)
            MPI_Comm_free(&in_group_comm);
    } else if (size <= HBCAST_MIN_PROCS || validate_input(num_groups, size)) {
        /*if (debug == 2)
           MPIX_Get_property(comm_world, MPIDO_RECT_COMM, &(bcast_response.rec_comm_world));
         */ 
        fprintf(stdout, "Using non-hierarchical bcast\n");
        hpnla_bcast(buffer, count, datatype, root, comm_world, alg);
    } else if (!validate_input(num_groups, size)) {
        /*TODO*/
        fprintf(stdout, "Wrong number of groups\n");
        bcast_response.err_flag = ERR_GROUPS;
    }

    return bcast_response;
}





/*
 * TODO:  It works only for root=0. Fix it to work with any possible root.
 */

/*t_bcast_response MPI_HBcast_Opt(void *buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm_world, MPI_Comm in_group_comm,
        MPI_Comm out_group_comm,  int num_groups, int rec, int alg, int debug) {
    int pg;
    int rank;
    int size;
    int stride;
    int root_inside;
    int root_outside;
    int my_group;

    t_bcast_response bcast_response = {
        .rec_out_group_comm = -1,
        .rec_in_group_comm = -1,
        .rec_comm_world = -1,
        .err_flag = 0
    };


    MPI_Comm_rank(comm_world, &rank);
    MPI_Comm_size(comm_world, &size);
    if (size == 1) return bcast_response;

    //TODO make num_groups configurable
    if (size > HBCAST_MIN_PROCS && validate_input(num_groups, size)) {

        //Start broadcast between groups
        if (out_group_comm != MPI_COMM_NULL && rec != -1) { // if rec == -1 then broadcast only inside
            int out_size = 0;
            int out_rank = -1;
            root_outside = root; //  root / pg; 
            hpnla_bcast(buffer, count, datatype, root_outside, out_group_comm, alg);          
        }

        root_inside = root; 
        
        //Start broadcast inside groups
        switch (rec) {
            case 3: // more levels of hierarchy, until size < HBCAST_MIN_PROCS
                MPI_HBcast_Opt(buffer, count, datatype, root_inside, in_group_comm, num_groups, 3, alg, debug);
                break;
            case 2: // two levels of hierarchy. Try making hierarchy further inside groups
                MPI_HBcast_opt(buffer, count, datatype, root_inside, in_group_comm, num_groups, 1, alg, debug);
                break;
            case 1: // 1 level of hierarchy
                hpnla_bcast(buffer, count, datatype, root_inside, in_group_comm, alg);
                break;
            case -1:
                // Just to see if broadcast inside groups is better than that of with MPI_COMM_WORLD
                hpnla_bcast(buffer, count, datatype, root_inside, in_group_comm, alg);
                break;
            default: // e.g. -2
                // Don't broadcast inside groups. Just to see if broadcast between groups is better than that of with MPI_COMM_WORLD
                break;
        }

    } else if (size <= HBCAST_MIN_PROCS || validate_input(num_groups, size)) {
        fprintf(stdout, "Using non-hierarchical bcast\n");
        hpnla_bcast(buffer, count, datatype, root, comm_world, alg);
    } else if (!validate_input(num_groups, size)) {
        fprintf(stdout, "Wrong number of groups\n");
        bcast_response.err_flag = ERR_GROUPS;
    }

    return bcast_response;
}

*/
