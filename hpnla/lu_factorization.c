/*!
 * Block-cylic distributed LU factorization 
 * algorithm with MPI
 *
 * Authors: Khalid Hasanov
 *
 */

#include "tools/hpnla_debug.h"
#include "tools/hpnla_timer.h"
#include "communication/hpnla_bcast.h"
#include "lu_factorization.h"
#include <stdlib.h>
#ifdef HDNLA_SMPI
#include <smpi.h>
#else
#define SMPI_SHARED_FREE free
#endif



double lu_factorize(LU_data* lu_data, Debug_data* debug_data) {
	double *a;
    int i, j;
    int ch;
    
    int me = lu_data->rank;
    int nproc = lu_data->nproc;
    int n = lu_data->n;
    int m = lu_data->m;
    int block_size = lu_data->block_size;
    int prow = lu_data->prow;
    int pcol = lu_data->prow;
    
    /*
     * Flags for different steps. Set zero to disable the corresponding step.
     */
    int with_partial_pivoting = lu_data->with_partial_pivoting;
    int print_proc_view = debug_data->print_proc_view;
    int do_ugly = debug_data->do_ugly_init;
    int print_initial_matrix = debug_data->print_initial_matrix;

    int step_pivoting = debug_data->step_pivoting;
    int step_dvide = debug_data->step_dvide;
    int step_l0_u1 = debug_data->step_l0_u1;
    int step_bcast_l_u = debug_data->step_bcast_l_u;
    int step_dgemm = debug_data->step_dgemm;
    int step_print_final = debug_data->step_print_final;
    int step_print_pvt = debug_data->step_print_pvt;

    

    if (me == 0) {
        printf("\n Block-Cyclic Distributed Dense LU Factorization\n");
        printf("     %d by %d Matrix\n", n, m);
        printf("     %d Processors\n", nproc);
        printf("     %d by %d Element Blocks\n", block_size, block_size);
        printf("\n");
    }


    {// Validation block 

        if (block_size < 2) {
            if (!me)
                fprintf(stderr, "Block size should be at least 2.\n");
            MPI_Finalize();
            exit(0);
        }


        if (n != m) {
            if (!me)
                fprintf(stderr, "Only square matrix is supported.\n");
            MPI_Finalize();
            exit(0);
        }

        if (n % block_size) {
            if (!me)
                fprintf(stderr, "n should be divisible by block_size.\n");
            MPI_Finalize();
            exit(0);
        }


        if (prow * pcol != nproc) {
            if (!me)
                fprintf(stderr, "num_rows * num_cols should be equal to nproc.\n");
            MPI_Finalize();
            exit(0);
        }

    } //End of the validation block


    int nblocks = n / block_size;

    int row = me / prow;
    int col = me % pcol;

    int n_loc = n / prow;
    int m_loc = m / pcol;

    MPI_Comm row_comm, col_comm;

    /* color by row, rank by column */
    MPI_Comm_split(lu_data->comm, row, col, &row_comm);
    /* color by column, rank by row */
    MPI_Comm_split(lu_data->comm, col, row, &col_comm);

    proc_grid diagonal_grid;
    

    if (print_proc_view)
        if (me == 0) {
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    diagonal_grid = bc_block_owner(i, j, block_size, block_size, prow, pcol);
                    fprintf(stdout, "P%d (%d, %d) ", diagonal_grid.pid, diagonal_grid.pi, diagonal_grid.pj);
                }
                fprintf(stdout, "\n");
            }
        }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);


    // Do some ugly hardcoded matrix init for p=4, n=8, m=8
    matrix_initialisation(&a, n_loc, m_loc, row, col, me, matrix_initialisation);

    if (print_initial_matrix) {
        print_matrix(0, 0, me, n_loc, m_loc, a, "a");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    int ib = 0;
    int r;

    struct {
        double val;
        int pvtline;
    } z, y;

    dbl_twoindex local, global;

    /* create our new data type */
    MPI_Datatype mpi_dbl_twoindex;
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_INT, MPI_INT};
    MPI_Aint disps[3] = {
        offsetof(dbl_twoindex, val),
        offsetof(dbl_twoindex, rank),
        offsetof(dbl_twoindex, posn)
    };
    int lens[3] = {1, 1, 1};
    MPI_Type_create_struct(3, lens, disps, types, &mpi_dbl_twoindex);
    MPI_Type_commit(&mpi_dbl_twoindex);

    /* create our operator */
    MPI_Op mpi_maxloc_dbl_twoindex;
    MPI_Op_create(maxloc_dbl_twoindex, 1, &mpi_maxloc_dbl_twoindex);


    //TODO: make pi, pj parametric by using offset
    int pi = 0;
    int pj = 0;

    int debug_rank = -2;

    if (me == debug_rank) {
        int loop = 1;
        char hostname[256];
        gethostname(hostname, sizeof (hostname));
        fprintf(stdout, "me=%d, PID = %d on %s ready for attach\n", me, getpid(), hostname);

        while (loop == 1);
    }

    int unprocessed_row = 0;
    int unprocessed_col = 0;

    int max_owner = -1;
    // int* pivot_rows = (int *) malloc(n * sizeof (int));
    // double* pivot_elms = (double *) malloc(n * sizeof (double));

    /*
     * The number of steps is n/b
     */
    for (i = 0; i < nblocks; i++) {
      //  if (i>0) continue;
        for (ib = 0; ib < block_size; ib++) {
            //if (i==0 && ib>0) continue;
            diagonal_grid = bc_block_owner(i * block_size + ib, i * block_size + ib,
                    block_size, block_size,
                    prow, pcol);

            if (col == pj) {
                if (with_partial_pivoting) {
                    if (unprocessed_row < n_loc) {
                        r = max_col_loc(a, n_loc, unprocessed_row + ib * (me == diagonal_grid.pid), unprocessed_col + ib);
                        if (r == -1) {
                            fprintf(stderr, "Something is wrong, me=%d, i=%d, ib=%d, unprocessed_col=%d, unprocessed_row=%d\n",
                                    me, i, ib, unprocessed_col, unprocessed_row);
                            exit(-1);
                        }

                        z.pvtline = r;
                        z.val = fabs(a[r * n_loc + unprocessed_col + ib]);
                    } else { //it should not happen
                        //   fprintf(stdout, "I've finished my part of the matrix! i = %d, ib = %d, me = %d\n", i, ib, me);
                        z.pvtline = r;
                        z.val = -86;
                        // This number will not be picked anyway, as we take absolute values in all other cases
                    }

                    local.val = z.val;
                    local.posn = r;
                    local.rank = me;

                    MPI_Allreduce(&local, &global, 1, mpi_dbl_twoindex, mpi_maxloc_dbl_twoindex, col_comm);

                    if (y.val == 0) {
                        if (unprocessed_row < n_loc && unprocessed_col < n_loc)
                            fprintf(stderr,
                                "I have a singular matrix.[i=%d, ib=%d, me=%d, y.pvtline=%d, y.val=%f, unprocessed_row=%d]\n",
                                i, ib, me, y.pvtline, y.val, unprocessed_row);

                        MPI_Finalize();
                        exit(-1);
                    }

                    max_owner = global.rank;
                    y.val = global.val;
                    y.pvtline = global.posn;

                    //TODO: don't do some ugly broadcast just to find out sign of the abs(max) element ??
                    int max_sign = -1;
                    if (me == max_owner && y.val == a[y.pvtline * n_loc + unprocessed_col + ib])
                        max_sign = 1;

                    MPI_Bcast(&max_sign, 1, MPI_INT, max_owner / prow, col_comm);
                    y.val *= max_sign;

                    MPI_Barrier(col_comm);
                    MPI_Barrier(col_comm);
                } else { // We just take the current diagonal element 
                    y.pvtline = unprocessed_row + ib * (me == diagonal_grid.pid);
                    y.val = a[y.pvtline * n_loc + unprocessed_col + ib];
                    MPI_Bcast(&y, 1, MPI_DOUBLE_INT, diagonal_grid.pi, col_comm);
                    fprintf(stdout, "Running without pivoting: [i:%d, ib:%d, me:%d, y.val:%f]\n", i, ib, me, y.val);
                }
            } // End of (col == pj)

            // Store global pivot row indexes
            if (me == diagonal_grid.pi && with_partial_pivoting) {
                //TODO: saving pivot info and printing it is not correct
                //int global_pivot_row = (max_owner / prow + prow * (y.pvtline / block_size) )* block_size + y.pvtline % block_size;
                int global_pivot_row = get_global_row(y.pvtline / block_size, y.pvtline % block_size, max_owner / prow, prow, block_size);
                if (step_print_pvt)
                    printf("i=%d, ib=%d, global_pivot_row=%d, y.val=%f, me=%d\n", i, ib, global_pivot_row, y.val, me);
            }


            if (step_pivoting) { // Block of broadcasting and exchanging pivot

                if (with_partial_pivoting) {
                    MPI_Bcast(&y, 1, MPI_DOUBLE_INT, diagonal_grid.pi, row_comm);
                    MPI_Bcast(&max_owner, 1, MPI_INT, diagonal_grid.pi, row_comm);

                    if (row == max_owner / prow && row == diagonal_grid.pi && y.pvtline == diagonal_grid.loc_i) {
                        // do nothing

                    } else if (row == max_owner / prow && row == diagonal_grid.pi) {
                        // local exchange                    
                        swap_with_pivot_row(diagonal_grid.loc_i, y.pvtline, n_loc, a);

                    } else if (row == max_owner / prow && row != diagonal_grid.pi) {
                        // exchange with the owner
                        MPI_Status status;

                        // I think I need to use allocate a buffer If I want to use a non-blocking send/receive 
                        MPI_Send(&a[ y.pvtline * n_loc], n_loc, MPI_DOUBLE, diagonal_grid.pi, PIVOT_TAG, col_comm);
                        MPI_Recv(&a[ y.pvtline * n_loc], n_loc, MPI_DOUBLE, diagonal_grid.pi, NON_PIVOT_TAG, col_comm, &status);

                    } else if (row != max_owner / prow && row == diagonal_grid.pi) {
                        /*
                         * Here we receive before we send our own data therefore we need to save it in a temporary buffer
                         */
                        double* current_row_buffer = (double *) malloc(n_loc * sizeof (double));
                        MPI_Status status;

                        MPI_Recv(current_row_buffer, n_loc, MPI_DOUBLE, max_owner / prow, PIVOT_TAG, col_comm, &status);
                        MPI_Send(&a[ diagonal_grid.loc_i * n_loc], n_loc, MPI_DOUBLE, max_owner / prow, NON_PIVOT_TAG, col_comm);

                        int il;
                        for (il = 0; il < n_loc; il++)
                            a[ diagonal_grid.loc_i * n_loc + il] = current_row_buffer[il];

                        free(current_row_buffer);

                    } else {
                        // printf("I am not interesting, leave me alone! "
                        //         "me=%d, max_owner=%d, grid.pid=%d, z.pvtline=%d, y.pvtline=%d, y.val=%f, i + ib=%zu\n",
                        //         me, max_owner, grid.pid, z.pvtline, y.pvtline, y.val, i + ib);                     
                    }

                }

                //Step 3: A(i+1:n, i) = A(i+1:n, i) / A(i, i), ie. calculation of L
                //TODO 3, move it into a function
                //Step 4: A(i+1:n, i+1:end) -= A(i+1:n, 1) * A(i, i+1:end)
                if (step_dvide) {
                    if (col == pj && (unprocessed_row < n_loc && unprocessed_col < n_loc)) {

                        double factor;
                        int row_to_dvide;
                        int col_to_factor;
                        int ipr;
                        int ifc;

                        row_to_dvide = (unprocessed_row + (me == diagonal_grid.pid)) + ib * (me == diagonal_grid.pid);

                        for (ipr = row_to_dvide; ipr < n_loc; ipr++) {
                            factor = a[ipr * n_loc + (unprocessed_col + ib)] / (double) y.val;
                            a[ipr * n_loc + (unprocessed_col + ib)] = factor;
                            col_to_factor = unprocessed_col + (col == diagonal_grid.pj) + ib;

                          //  printf("i=%d, ib=%d, me=%d, col=%d, col_to_factor=%d\n", i, ib, me, (unprocessed_col + ib), col_to_factor);

                            //Step 4: There is no step 4 in lawn58_lapack ??
                            for ( ifc = col_to_factor; ifc < unprocessed_col + block_size; ifc++ )
                                 a[ipr * n_loc + ifc] -= factor * a[unprocessed_row + ifc]; 
                             
                        }
                    }
                } // End of divide blocks

            } // End of pivot and bcast block
        }


        if (step_l0_u1) { // Start block of calculation of L0 and U1

            if (row == pi) {
                /* 
                 * TODO: make it general.
                 * Calculation of max number of non-zero elements in LL
                 * assumes the matrix A and blocks are square. Then the 
                 * number of elements in the diagonal is subtracted.
                 */
                int max_num_nnz = (block_size) * (block_size + 1) / 2 - block_size;
                double *l0_partial_buffer = (double*) malloc(max_num_nnz * sizeof (double));

                int l0_i, l0_j, l0_k = 0;
                if (col == pj) {
                    for (l0_i = unprocessed_row; l0_i < unprocessed_row + block_size; l0_i++) {
                        for (l0_j = unprocessed_col; l0_j < unprocessed_col + block_size; l0_j++) {
                            if (l0_i > l0_j) {
                                l0_partial_buffer[l0_k++] = a[l0_i * n_loc + l0_j];
                            }
                        }
                    }
                }

                int u1_start_col = (me == diagonal_grid.pid) ? unprocessed_col + block_size : unprocessed_col;
                double* u1_buffer = (double*) malloc(block_size * (n_loc - u1_start_col) * sizeof (double));


                //TODO. should I check col==pj ??
                int u1_i, u1_j, u1_k = 0;
                for (u1_i = unprocessed_row; (u1_i < unprocessed_row + block_size); u1_i++) {
                    for (u1_j = u1_start_col; u1_j < n_loc; u1_j++) {
                        u1_buffer[u1_k] = a[u1_i * n_loc + u1_j];
                        u1_k++;
                    }
                }


                //TODO:  start benchmarking bcast

                /*
                 * Bcast L0 right
                 */
                MPI_Bcast(l0_partial_buffer, max_num_nnz, MPI_DOUBLE, diagonal_grid.pi, row_comm);

                //TODO:  end benchmarking bcast



                double *l0_buffer = (double*) malloc(block_size * block_size * sizeof (double));

                l0_k = 0;
                for (l0_i = 0; l0_i < block_size; l0_i++) {
                    for (l0_j = 0; l0_j < block_size; l0_j++) {
                        if (l0_i == l0_j)
                            l0_buffer[l0_i * block_size + l0_j] = 1.0;
                        else if (l0_i < l0_j)
                            l0_buffer[l0_i * block_size + l0_j] = 0.0;
                        else {
                            l0_buffer[l0_i * block_size + l0_j] = l0_partial_buffer[l0_k];
                            l0_k++;
                        }
                    }
                }

                free(l0_partial_buffer);

                int order = CblasRowMajor;
                int side = CblasLeft;
                int uplo = CblasLower;
                int trans = CblasNoTrans;
                int diag = CblasUnit;
                int alpha = 1.0;
                int lda = block_size;
                int ldb = n_loc - u1_start_col;


                //TODO. ldb==0 means that all the columns are processed already
                if (ldb != 0) {
                    //	printf("BEFOR, i=%d, me=%d, u1[0]=%f, u1[1]=%f, u1[2]=%f, u1[3]=%f\n", i, me, u1_buffer[0], u1_buffer[1], u1_buffer[2], u1_buffer[3]);
                    //block_size * (n_loc - u1_start_col)					
                    // L0xU1 = C
                    cblas_dtrsm(order, side, uplo, trans, diag, lda, ldb, alpha, l0_buffer, lda, u1_buffer, ldb);

                    u1_k = 0;
                    for (u1_i = unprocessed_row; (u1_i < unprocessed_row + block_size); u1_i++) {
                        for (u1_j = u1_start_col; u1_j < n_loc; u1_j++) {
                            a[u1_i * n_loc + u1_j] = u1_buffer[u1_k];
                            u1_k++;
                        }
                    }

                    //printf("AFTER, i=%d, me=%d, u1[0]=%f, u1[1]=%f, u1[2]=%f, u1[3]=%f\n", i, me, u1_buffer[0], u1_buffer[1], u1_buffer[2], u1_buffer[3]);     
                }

                free(l0_buffer);
                free(u1_buffer);

            }

        } // end of block step_l0_u1


        /*
         *  Start horizontal and vertical bcasts then dgemm
         */
        if (step_bcast_l_u) {
            int u_start_col = (col == diagonal_grid.pj) ? unprocessed_col + block_size : unprocessed_col;
            double* u_buffer = (double*) malloc(block_size * (n_loc - u_start_col) * sizeof (double));

            /*
             * Bcast A(ib:end, end+1:n) down 
             */
            if (u_start_col < n_loc) { //If not, it means there is no one to bcast ??                
                int u_i, u_j, u_k = 0;
                for (u_i = unprocessed_row; (u_i < unprocessed_row + block_size); u_i++) {
                    for (u_j = u_start_col; u_j < n_loc; u_j++) {
                        u_buffer[u_k] = a[u_i * n_loc + u_j];
                        u_k++;
                    }
                }

                //TODO:  start benchmarking bcast 
                MPI_Bcast(u_buffer, block_size * (n_loc - u_start_col), MPI_DOUBLE, diagonal_grid.pi, col_comm);
                //TODO:  end benchmarking bcast
            }

            int l_start_row = (row == diagonal_grid.pi) ? unprocessed_row + block_size : unprocessed_row;
            double* l_buffer = (double*) malloc((n_loc - l_start_row) * block_size * sizeof (double));


            /*
             * Bcast A(end+1:n, ib:end) right
             */
            if (l_start_row < n_loc) {
                int l_i, l_j, l_k = 0;
                for (l_i = l_start_row; l_i < n_loc; l_i++) {
                    for (l_j = unprocessed_col; (l_j < unprocessed_col + block_size); l_j++) {
                        l_buffer[l_k] = a[l_i * n_loc + l_j];
                        l_k++;
                    }
                }

                //TODO:  start benchmarking bcast
                MPI_Bcast(l_buffer, (n_loc - l_start_row) * block_size, MPI_DOUBLE, diagonal_grid.pj, row_comm);
                //TODO:  end benchmarking bcast    
            }


            /*
             * Eliminate A(end+1:n, end+1:n)
             */
            if (step_dgemm) {
                if (l_start_row < n_loc && u_start_col < n_loc) {
                    int alpha = -1;
                    int beta = 1;

                    // E = E - L1 x U1
                    /* C = alpha x A x B + beta x C
                     * C:  m x n
                     * A:  m x k
                     * B:  k x n
                     */
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            // m, n, k
                            n_loc - l_start_row, n_loc - u_start_col, block_size,
                            alpha,
                            // A
                            l_buffer,
                            // lda
                            block_size,
                            // B
                            u_buffer,
                            // ldb
                            n_loc - u_start_col,
                            beta,
                            // C
                            &a[ l_start_row * n_loc + u_start_col ],
                            // ldc);
                            n_loc);

                }

            } // End of block step_dgemm

            free(u_buffer);
            free(l_buffer);

        } // End of block step_bcast_l_u


        if (col == pj) {
            unprocessed_col += block_size;
        }


        if (row == pi) {
            unprocessed_row += block_size;
        }

        pj = (pj + 1) % pcol;
        pi = (pi + 1) % prow;

    }

    if (step_print_final) {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if (!me)
            printf("At the end\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        print_matrix(0, 0, me, n_loc, m_loc, a, "a");

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
    } // End of block step_print_final

    // free(pivot_rows);
    free(a);

    MPI_Op_free(&mpi_maxloc_dbl_twoindex);
    MPI_Type_free(&mpi_dbl_twoindex);

    MPI_Finalize();

    return 0;
	
}









int get_global_row(int l, int x, int rank, int prow, int brow) {
    return (l * prow + rank) * brow + x;
}

int get_global_col(int m, int y, int rank, int pcol, int bcol) {
    return (m * pcol + rank) * bcol + y;
}

/* TODO:
 * For now, brow=bcol. Make nb and mb separately
 */
proc_grid bc_block_owner(int i, int j, int brow, int bcol, int prow, int pcol) {
    proc_grid grid;

    grid.pi = (int) (floor(i / brow)) % prow;
    grid.pj = (int) (floor(j / bcol)) % pcol;
    grid.pid = grid.pi * pcol + grid.pj;

    grid.l = (int) (floor(i / (prow * brow)));
    grid.m = (int) (floor(j / (pcol * bcol)));

    grid.x = (int) (i % brow);
    grid.y = (int) (j % bcol);

    grid.loc_i = (int) (grid.l * brow + grid.x);
    grid.loc_j = (int) (grid.m * bcol + grid.y);

    //TODO
    return grid;
}

/* TODO */
void swap_with_pivot_row(int source_row, int pivot_row, int col_size, double* a) {
    int i;
    double* temp_row = malloc(col_size * sizeof (double));

    if (source_row != pivot_row) {
        for (i = 0; i < col_size; i++) {
            temp_row[i] = a[source_row * col_size + i];
            a[source_row * col_size + i] = a[pivot_row * col_size + i];
            a[pivot_row * col_size + i] = temp_row[i];
        }

        free(temp_row);
    }

}

void print_matrix(int is, int ib, int rank, int n, int m, double* a, char* matrix_name) {
    int i, j;

    for (i = 0; i < n; i++) {
        printf("me=%d -> ", rank);
        for (j = 0; j < m; j++) {
            printf("%s[%d][%d]=%f ", matrix_name, i, j, a[i * n + j]);
        }
        printf("\n");
    }

}

void matrix_initialisation(double** p_a, int n_loc, int m_loc, int row, int col, int rank, int do_ugly) {
    int x, z;
    int lda = m_loc;
    double *a;


    if (do_ugly) {
		 if (me==0) {
			   a[0]=6.0; a[1]=9.0; a[2]=9.0; a[3]=3.0;
			   a[4]=2.0; a[5]=4.0; a[6]=2.0; a[7]=2.0;
			   a[8]=1.0; a[9]=9.0; a[10]=7.0; a[11]=3.0;
			   a[12]=2.0; a[13]=7.0; a[14]=0.0; a[15]=9.0;
	     } else if (me==1) {
			   a[0]=6.0; a[1]=7.0; a[2]=4.0; a[3]=5.0;
			   a[4]=5.0; a[5]=6.0; a[6]=8.0; a[7]=6.0;
			   a[8]=6.0; a[9]=6.0; a[10]=3.0; a[11]=6.0;
			   a[12]=2.0; a[13]=10.0; a[14]=6.0; a[15]=5.0;		   
	     } else if (me==2) {
			   a[0]=2.0; a[1]=1.0; a[2]=11.0; a[3]=9.0;
			   a[4]=0.0; a[5]=10.0; a[6]=8.0; a[7]=6.0;
			   a[8]=9.0; a[9]=1.0; a[10]=5.0; a[11]=3.0;
			   a[12]=4.0; a[13]=1.0; a[14]=8.0; a[15]=1.0;
	     } else if (me==3) {
			   a[0]=4.0; a[1]=7.0; a[2]=12.0; a[3]=5.0;
			   a[4]=2.0; a[5]=11.0; a[6]=3.0; a[7]=6.0;
			   a[8]=11.0; a[9]=6.0; a[10]=9.0; a[11]=9.0;
			   a[12]=1.0; a[13]=6.0; a[14]=8.0; a[15]=12.0;
	     }	   	   	   
	} else {
		srand(time(NULL));

		if (p_a != NULL) {
			a = malloc(sizeof (double) * n_loc * m_loc);
			if (a == 0) {
				perror("Error allocation Matrix A");
				exit(-1);
			}

			*p_a = a;

			//initialization of the matrices
			for (x = 0; x < n_loc; x++) {
				for (z = 0; z < m_loc; z++) {
					a[x * lda + z] = (double) (rand() % (rank + 10)); //(double) (x + col * n_loc + z);
				}
			}
		}
    }

}

int max_col_loc(double *a, int num_rows, int row, int col) {
    int i;
    int r = row;

    double max = fabs(a[row * num_rows + col]);
    for (i = row + 1; i < num_rows; i++) {
        if (fabs(a[i * num_rows + col]) > max) {
            max = a[i * num_rows + col];
            r = i;
        }
    }

    return r;
}

/*
 * This custom MPI_Op funciton was adapted from 
 *   http://stackoverflow.com/questions/9285442/mpi-get-processor-with-minimum-value
 */
void maxloc_dbl_twoindex(void *in, void *inout, int *len, MPI_Datatype *type) {
    /* ignore type, just trust that it's our dbl_twoindex type */
    dbl_twoindex *invals = in;
    dbl_twoindex *inoutvals = inout;

    int i;

    for (i = 0; i<*len; i++) {
        if (invals[i].val >= inoutvals[i].val) {
            inoutvals[i].val = invals[i].val;
            inoutvals[i].rank = invals[i].rank;
            inoutvals[i].posn = invals[i].posn;
        }
    }

    return;
}

/*
 * TODO: use daxpy from BLAS
 */
void daxpy(double *a, double *b, int n, double alpha) {
    int i;

    for (i = 0; i < n; i++) {
        a[i] += alpha * b[i];
    }
}





