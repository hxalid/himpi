/* 
 * File:   lu_factorize.h
 * Author: khalid
 *
 * Created on 25 November 2013, 14:29
 */

#ifndef LU_FACTORIZE_H
#define	LU_FACTORIZE_H

#include "cblas_wrappers/hpnla_cblas.h"

#ifdef	__cplusplus
extern "C" {
#endif

    
     typedef struct LU_data {
        hpnla_gemm* gemm;
        size_t m_global;
        size_t k_global;
        size_t n_global;
        size_t m;
        size_t n;
        size_t k_a;
        size_t k_b;
        size_t Block_size;
        size_t key;
        size_t row;
        size_t col;
        int bcast_algorithm;
        int distribution;       // not used for now

    } LU_data;

    typedef struct Platform_data {
        size_t size_row;
        size_t size_col;
        int my_rank;
        int nb_proc;
        size_t nb_requested_proc;
        MPI_Comm comm;
    } Platform_data;


    double lu_factorize(LU_data* lu_data, Platform_data* platform_data);

    int validate_input(LU_data* lu_data, Platform_data* platform_data);


    



#ifdef	__cplusplus
}
#endif

#endif	/* LU_FACTORIZE_H */

