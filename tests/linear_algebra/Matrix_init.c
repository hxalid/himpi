#include "config.h"
#include "Matrix_init.h"
#include "tools/hpnla_debug.h"
#include <math.h>
#define _unused(x) ((void)x)
#ifdef HDNLA_SMPI
#include <smpi.h>
#else
#define SMPI_SHARED_MALLOC malloc
#endif

void matrices_initialisation(double ** p_a, double ** p_b, double ** p_c,
        size_t m, size_t k_a, size_t k_b, size_t n,
        size_t row, size_t col) {

    size_t x, y, z;
    size_t lda = k_a;
    size_t ldb = n;
    size_t ldc = n;
    double *a, *b, *c;
    _unused(row);

    //TODO explain the parameter -1
    a = SMPI_SHARED_MALLOC(sizeof (double) * m * k_a);


    if (a == 0) {
        perror("Error allocation Matrix A");
        exit(-1);
    }

    //TODO explain the parameter -1
    b = (double *) SMPI_SHARED_MALLOC(sizeof (double) * k_b * n
            );

    if (b == 0) {
        perror("Error allocation Matrix B");
        exit(-1);
    }

    //TODO explain the parameter 
    c = (double *) SMPI_SHARED_MALLOC(sizeof (double) * m * n
            );
    if (c == 0) {
        perror("Error allocation Matrix C");
        exit(-1);
    }

    *p_a = a;
    *p_b = b;
    *p_c = c;
    
    
    //initialization of the matrices
    for (x = 0; x < m; x++) {
        for (z = 0; z < k_a; z++) {
#ifdef SIMPLE_MATRIX
            a[x * lda + z] = 1;
#else
            a[x * lda + z] = (double) (z + col * n);
  //          printf("a[%d][%d]=%f\t", x, z, a[x * lda + z]);
#endif
        }
 //       printf("\n");
    }
    
 //   printf("*********m=%d, n=%d, k_b=%d, k_a=%d********\n",m, n, k_b, k_a);
    for (z = 0; z < k_b; z++) {
        for (y = 0; y < n; y++) {
#ifdef SIMPLE_MATRIX
            b[z * ldb + y] = 1;
#else
            b[z * ldb + y] = (double) (y);
        //    printf("b[%d][%d]=%f\t", z, y, b[z * ldb + y]);
#endif
        }
    //    printf("\n");
    }
    
    for (x = 0; x < m; x++) {
        for (y = 0; y < n; y++) {
            c[x * ldc + y] = 0;
        }
    }
}

void matrices_allocation(double ** p_a, double ** p_b, double ** p_c,
        size_t m, size_t k_a, size_t k_b, size_t n) {

    double * a, *b, *c;

    //TODO explain the parameter -1
    a = (double *) SMPI_SHARED_MALLOC(sizeof (double) * m * k_a
            );

    if (a == 0) {
        perror("Error allocation Matrix A");
        exit(-1);
    }

    //TODO explain the parameter -1
    b = (double *) SMPI_SHARED_MALLOC(sizeof (double) * k_b * n
            );

    if (b == 0) {
        perror("Error allocation Matrix B");
        exit(-1);
    }

    //TODO explain the parameter 
    c = (double *) SMPI_SHARED_MALLOC(sizeof (double) * m * n
            );
    if (c == 0) {
        perror("Error allocation Matrix C");
        exit(-1);
    }

    *p_a = a;
    *p_b = b;
    *p_c = c;

}

void blocks_initialisation(double ** p_a_local, double ** p_b_local,
        size_t m, size_t B_k, size_t n) {
    size_t x, y, z;
    size_t lda = B_k;
    size_t ldb = n;
    double * a_local, *b_local;

    //TODO explain the parameter -1
    a_local = (double *) SMPI_SHARED_MALLOC(sizeof (double) * m * B_k
            );

    if (a_local == 0) {
        perror("Error allocation Matrix A");
        exit(-1);
    }

    //TODO explain the parameter -1
    b_local = (double *) SMPI_SHARED_MALLOC(sizeof (double) * B_k * n
            );

    if (b_local == 0) {
        perror("Error allocation Matrix B");
        exit(-1);
    }

    *p_a_local = a_local;
    *p_b_local = b_local;

    // initialization of the matrices
    for (x = 0; x < m; x++) {
        for (z = 0; z < B_k; z++) {
            a_local[x * lda + z] = 0.0;
        }
    }
    for (z = 0; z < B_k; z++) {
        for (y = 0; y < n; y++) {
            b_local[z * ldb + y] = 0.0;
        }
    }
}

void check_reception_a(double *a_local,
        size_t m, size_t B_k, size_t row, size_t col) {
    size_t x, z;
    size_t lda = B_k;
    for (x = 0; x < m; x++) {
        for (z = 0; z < B_k; z++) {
            if (a_local[x * lda + z] != 2.0) {
                debug_print(2, row, col,
                        "the a_local is not good at all (%zu,%zu): %f\n",
                        x, z, a_local[x * lda + z]);
            }
        }
    }
}

void check_reception_b(double *b_local,
        size_t B_k, size_t n, size_t row, size_t col) {
    size_t y, z;
    size_t ldb = n;
    for (z = 0; z < B_k; z++) {
        for (y = 0; y < n; y++) {
            if (b_local[z * ldb + y] != 1.0) {
                debug_print(0, row, col,
                        "the b_local is not good at all (%zu,%zu): %f\n",
                        z, y, b_local[z * ldb + y]);
                break;
            }
        }
        if (b_local[z * ldb + y] != 1.0) {
            break;
        }
    }
}

void check_result(double *c, double *a, double *b,
        size_t m, size_t n, size_t k_a, size_t k_b,
        size_t row, size_t col,
        size_t size_row, size_t size_col) {
    size_t x, y;
    size_t ldc = n;
    _unused(a);
    _unused(b);
    _unused(k_b);
    _unused(k_a);
    _unused(row);
    _unused(col);
    _unused(size_row);
    /* these variable could be use to check the result in function of the
     * matrix initialization */


    /*Display for checking */
#ifdef SIMPLE_MATRIX
    printf("Value get : %f excepted %zu multiply by y\n", c[((int) m / 2) * ldc + 1], size_row * k_a);
#else
    printf("Value get : %f excepted %zu multiply by y\n", c[((int) m / 2) * ldc + 1], 1 * (size_col * m)*((size_col * m) - 1) / 2);
#endif
    for (x = 0; x < m; x++) {
        for (y = 0; y < n; y++) {
            /* WARNING this could be lead to some errors ( precision with double )*/
#ifdef SIMPLE_MATRIX
            if (fabs(c[x * ldc + y] - size_row * k_a) > 0.0000001)
#else
            if (fabs(c[x * ldc + y] - y * (size_col * m)*((size_col * m) - 1) / 2) > 0.0000001)
#endif
            {
#ifdef SIMPLE_MATRIX
                printf("%f\t%zu, y : %zu x : %zu \n",
                        c[x * ldc + y], size_row*k_a, y, x);
#else
                printf("%f\t%zu, y : %zu x : %zu \n",
                        c[x * ldc + y], y * (size_col * m)*((size_col * m) - 1) / 2, y, x);
#endif
                goto error_exit;
            }
        }
    }
    printf("result check: ok\n");
    return;
error_exit:
    printf("result check not ok\n"
            "WARNING the test could be lead to some "
            "errors ( precision with double )\n");
    return;
}


