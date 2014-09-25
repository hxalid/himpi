/* 
 * File:   lu_factorization.h
 * Author: Khalid Hasanov
 *
 * Created on 22 November 2013, 12:44
 */

#ifndef LU_FACTORIZATION_H
#define	LU_HFACTORIZATION_H

#include "cblas_wrappers/hpnla_cblas.h"

#define max(_a, _b) ( (_a) < (_b) ? (_b) : (_a) )
#define min(_a, _b) ( (_a) > (_b) ? (_b) : (_a) )

#define PIVOT_TAG 86
#define NON_PIVOT_TAG 88


#ifdef	__cplusplus
extern "C" {
#endif

    typedef struct LU_data {
		//
		int n;
		int m;
		int k;
		int block_size;
		//
		int nproc;
		int prow;
		int pcol;
		int rank;
		MPI_Comm comm;
		//
		int with_partial_pivoting;
		int bcast_algorithm;
		//
		hpnla_gemm* gemm;
	} LU_data;
	
	typedef struct Debug_data {		
		int print_proc_view;
		int print_initial_matrix;
		int do_ugly_init;		
		int step_pivoting;
		int step_dvide;
		int step_l0_u1;
		int step_bcast_l_u;
		int step_dgemm;
		int step_print_final;
		int step_print_pvt;
	} Debug_data;

	typedef struct proc_grid {
		int pi;
		int pj;
		int pid;
		int l;
		int m;
		int x;
		int y;
		int loc_i;
		int loc_j;
	} proc_grid;
	

	typedef struct dbl_twoindex_struct {
		double val;
		int rank;
		int posn;
	} dbl_twoindex;
	

	int get_global_row(int l, int x, int rank, int prow, int brow);

	int get_global_col(int m, int y, int rank, int pcol, int bcol);
	
	proc_grid bc_block_owner(int i, int j, int brow, int bcol, int prow, int pcol);
	
	void swap_with_pivot_row(int source_row, int pivot_row, int col_size, double* a);
  
    void print_matrix(int is, int ib, int rank, int n, int m, double* a, char* matrix_name);
    
    void matrix_initialisation(double** p_a, int n_loc, int m_loc, int row, int col, int rank);
    
    int max_col_loc(double *a, int num_rows, int row, int col);
    
    void maxloc_dbl_twoindex(void *in, void *inout, int *len, MPI_Datatype *type);
    

    double lu_factorize(LU_data* lu_data, Debug_data* debug_data);

    int validate_input(LU_data* lu_data, Debug_data* debug_data);




#ifdef	__cplusplus
}
#endif

#endif	/* LU_FACTORIZATION_H */

