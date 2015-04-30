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


//const char *HMPI_CONF_FILE_NAME = "fayil.conf";


int HMPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm, int rec, int alg) {
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


   // int num_groups = hmpi_get_num_groups(comm, HMPI_CONF_FILE_NAME);

    /*TODO make num_groups configurable*/
    if (comm_size > HBCAST_MIN_PROCS && validate_groups(num_groups, comm_size)) {
        pg = comm_size / num_groups;
        my_group = rank / pg;
        stride = root / pg;

        MPI_Comm_split(comm, (rank - my_group * pg == 0) ? 0 : MPI_UNDEFINED, rank, &out_group_comm);
        //  MPI_Comm_split(comm_world, (rank - my_group * pg == stride) ? 0 : MPI_UNDEFINED, rank, &out_group_comm);

#if(2 == DEBUG)
        MPIX_Get_property(in_group_comm, MPIDO_RECT_COMM, &(bcast_response.rec_in_group_comm));
#endif

        /*
         * Start broadcast between groups
         */
        if (out_group_comm != MPI_COMM_NULL && rec != -1) { // if rec == -1 then broadcast only inside
            int out_size = 0;
            int out_rank = -1;
            root_outside = root; //  root / pg; 
            hpnla_bcast(buffer, count, datatype, root_outside, out_group_comm, alg);

#if(2 == DEBUG)
            MPIX_Get_property(out_group_comm, MPIDO_RECT_COMM, &(bcast_response.rec_out_group_comm));
#endif

        }

        /*
         * Start broadcast inside groups
         */
        MPI_Comm_split(comm, my_group, rank, &in_group_comm);
        root_inside = root; // root_inside = stride;
        switch (rec) {
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
    } else if (comm_size <= HBCAST_MIN_PROCS || validate_groups(num_groups, comm_size)) {

#if(2 == DEBUG)
    	MPIX_Get_property(comm, MPIDO_RECT_COMM, &(bcast_response.rec_comm_world));
#endif
        
        fprintf(stdout, "Using non-hierarchical bcast\n");
        hpnla_bcast(buffer, count, datatype, root, comm, alg);
    } else if (!validate_groups(num_groups, comm_size)) {
        /*TODO*/
        fprintf(stdout, "Wrong number of groups\n");
    }

    return MPI_SUCCESS;
        
}




