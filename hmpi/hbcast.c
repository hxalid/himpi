#include <stdio.h>
#include <stdlib.h>

#include "hmpi.h"
#include "tools/utils.h"
#include "communication/hpnla_bcast.h"

#ifdef MPIX_H
#include <mpix.h>
#else 
#define MPIX_Get_property do{}while(0);
#define MPIDO_RECT_COMM 0
#endif


int hierarchical_broadcast(void *buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm, int num_groups, int num_levels, int alg_in, int alg_out) {
    int pg;
    int rank;
    int comm_size;
    int stride;
    int root_inside;
    int root_outside;
    int my_group;

    MPI_Comm in_group_comm, out_group_comm;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);
    
    if (comm_size == 1) return MPI_SUCCESS;

    if (comm_size > HBCAST_MIN_PROCS && validate_groups(num_groups, comm_size)) {
        pg = comm_size / num_groups;
        my_group = rank / pg;
        stride = root / pg;

        MPI_Comm_split(comm, (rank - my_group * pg == 0) ? 0 : MPI_UNDEFINED, rank, &out_group_comm);

#if(2 == DEBUG)
        MPIX_Get_property(in_group_comm, MPIDO_RECT_COMM, &(bcast_response.rec_in_group_comm));
#endif

        /*
         * Start broadcasting between groups.
         * If num_levels == -1 then broadcast only inside (just for debugging),
         * else if num_levels is 1 then do hierarchical bcast,
         * else broadcast only between groups (just for debugging).
         */
        if (out_group_comm != MPI_COMM_NULL && num_levels != -1) {
            int out_size = 0;
            int out_rank = -1;
            root_outside = root;
            hpnla_bcast(buffer, count, datatype, root_outside, out_group_comm, alg_out);

#if(2 == DEBUG)
            MPIX_Get_property(out_group_comm, MPIDO_RECT_COMM, &(bcast_response.rec_out_group_comm));
#endif

        }

        /*
         * Start broadcasting inside groups
         */
        MPI_Comm_split(comm, my_group, rank, &in_group_comm);
        root_inside = root;
        switch (num_levels) {
            case 1: // 1 level of hierarchy
                hpnla_bcast(buffer, count, datatype, root_inside, in_group_comm, alg_in);
                break;
            case -1:
                // Just to see if broadcast inside groups is better than that of with MPI_COMM_WORLD
                hpnla_bcast(buffer, count, datatype, root_inside, in_group_comm, alg_in);
                break;
            default: // e.g. -2
                // Don't broadcast inside groups. Just to see if broadcast between groups is better than that of with MPI_COMM_WORLD
                break;
        }

        if (out_group_comm != MPI_COMM_NULL)
            MPI_Comm_free(&out_group_comm);
        if (in_group_comm != MPI_COMM_NULL)
            MPI_Comm_free(&in_group_comm);
    } else if (comm_size <= HBCAST_MIN_PROCS || validate_groups(num_groups, comm_size)) {

#if(2 == DEBUG)
    	MPIX_Get_property(comm, MPIDO_RECT_COMM, &(bcast_response.rec_comm_world));
#endif
        
       // fprintf(stdout, "Using non-hierarchical bcast\n");
    	//TODO: alg_in is default is there is no hierarchy
        hpnla_bcast(buffer, count, datatype, root, comm, alg_in);
    } else if (!validate_groups(num_groups, comm_size)) {
        /*TODO*/
    	hpnla_bcast(buffer, count, datatype, root, comm, alg_in);
        fprintf(stdout, "Wrong number of groups: switching back to MPI_bcast\n");
    }

    return MPI_SUCCESS;
}

int HMPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm) {
	hmpi_conf my_conf = hmpi_get_my_conf(comm, HMPI_CONF_FILE_NAME);

	return hierarchical_broadcast(buffer, count, datatype, root, comm, my_conf.num_groups,
			my_conf.num_levels, my_conf.alg_in, my_conf.alg_out);
}


