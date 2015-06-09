/*!
 * 2.5D Block Matrix Multiplication example
 *
 * Authors: Quintin Jean-Noël
 */

#include "tools/hpnla_debug.h"
#include "tools/hpnla_timer.h"
#include "hmpi/mpi_bcast_algs.h"
#include "Matrix_init.h"
#include "Summa.h"
#include <stdlib.h>
#ifdef HDNLA_SMPI
#include <smpi.h>
#else
#define SMPI_SHARED_MALLOC malloc
#define SMPI_SHARED_FREE free
#endif


double two_dot_five(hpnla_gemm* gemm,
                    size_t m, size_t k, size_t n,
                    size_t Block_size, size_t group, size_t key,
                    size_t size_row, size_t size_col, size_t NB_groups ,
                    MPI_Comm world, int bcast_algorithm, int distribution){
  double *a,  *b,  *c, *res;
  /* Split the communicator into groups */

  /* Find out my identity in the default communicator */
  int myrank;
  int n_global=n;
//  int m_global=m;
//  int k_global=k;
  int NB_proc;
  int err;
  int useless = 0;

  init_timer();

  double time, communication_time = 0, summa_communication_t=0, summa_computation_t=0;
  struct timespec start_time, end_time; //time mesure
  struct timespec end_time_intern; //time mesure
  struct timespec start_time_reduce, end_time_reduce; //time mesure

  MPI_Comm my_world;

  if ( group >= NB_groups ){
    debug_print(1, (size_t) 0, (size_t) 0,
                "Not enough group NB_groups : %zu my group id : %zu\n",
                NB_groups, group);
    MPI_Comm_split(world, 0, key, &my_world);
    return -1;
  }else{
    MPI_Comm_split(world, 1, key, &my_world);
  }

  MPI_Comm_size (my_world, &NB_proc);

  if ( NB_proc < (int)(size_row*size_col*NB_groups) ){
    info_print(0, (size_t) 0, (size_t) 0,
                "Not enough processors NB_proc : %d required : %zu\n",
                NB_proc, size_row*size_col*NB_groups);
    return -1;
  }

  MPI_Comm group_comm;
  MPI_Comm_split(my_world, group, key, &group_comm);

  MPI_Comm_rank(group_comm, &myrank);
  MPI_Comm_size (group_comm, &NB_proc);
  /* for each group start the execution of his */

  NB_proc=size_row*size_col;
  size_t row = myrank / size_row;
  size_t col = myrank % size_row;


  /*-------------------------Check some mandatory conditions------------------*/
  size_t NB_Block = k / Block_size;
  if ( k % Block_size != 0 ){
    info_print(0, row, col,
                "The matrix size has to be proportionnal to the number\
                of blocks: %zu\n", NB_Block);
    return -1;
  }

  if ( size_row > NB_Block || size_col > NB_Block ){
    info_print(0, row, col,
                "Number of blocks is too small compare to the number of"
                " processors (%zu,%zu) in a row or a col (%zu)\n",
                size_col, size_row, NB_Block);
    return -1;
  }

  if ( NB_Block % size_row != 0 || NB_Block % size_col != 0){
    info_print(0, row, col,"The number of Block by processor is not an %s",
                "integer\n");
    return -1;
  }

  if(row >= size_col || col >= size_row){
    info_print(0, row, col, "I'm useless bye!!! col: %zu row: %zu, "
                "size_col: %zu , size_row: %zu \n",
                col,row,size_col,size_row);
    useless = 1;
  }


  if(useless == 1){
    /*----------------------Prepare the Communication Layer-------------------*/
    /* add useless processor on a new color to execute the matrix
     * multiplication with the other processors*/

    /* Split comm size_to row and column comms */
    MPI_Comm row_comm, col_comm, group_line;
    MPI_Comm_split(my_world, myrank, MPI_UNDEFINED, &group_line);
    /* color by row, rank by column */
    MPI_Comm_split(group_comm, size_row, MPI_UNDEFINED, &row_comm);
    /* color by column, rank by row */
    MPI_Comm_split(group_comm, size_col, MPI_UNDEFINED, &col_comm);
    /*------------------------Communication Layer can be used-----------------*/

    return 0;
  }

  debug_print(1, row, col, "I'm initialized col: %zu row: %zu, "
              "size_col: %zu , size_row: %zu, my rank: %d \n"
              ,col,row,size_col,size_row, myrank);



  /*------------------------Initialize the matrices---------------------------*/

  /* think about a common interface
   *  int pdgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans,
   *                 m, n, k, alpha, a, ia, ja, lda, b, ib, jb, ldb,
   *                 beta, c, ldc, Comm, rank );
   */

  /*------------------------Prepare the Communication Layer-------------------*/
  /* Split comm size_to row and column comms */
  MPI_Comm row_comm, col_comm, group_line;
  MPI_Comm_split(my_world, myrank, group, &group_line);
  /* color by row, rank by column */
  MPI_Comm_split(group_comm, row, col, &row_comm);
  /* color by column, rank by row */
  MPI_Comm_split(group_comm, col, row, &col_comm);
  /*-------------------------Communication Layer can be used------------------*/

  // matrix sizes
  m   = m   / size_col;
  n   = n   / size_row;
  size_t k_a = k / size_row;
  size_t k_b = k / size_col;


  /*only on the group 0*/
  if( group == 0 ) {
    matrices_initialisation(&a, &b, &c, m, k_a, k_b, n, row, col);
    if( NB_groups > 1 ) res = SMPI_SHARED_MALLOC( sizeof(double)*m*n );
  } else matrices_allocation(&a, &b, &c, m, k_a, k_b, n);


  /*-------------------Configuration for Summa algorihtm--------------------*/
  /*--------------------Allocation of matrices block-------------------------*/
  double *a_Summa, *b_Summa;
  blocks_initialisation(&a_Summa, &b_Summa, m, Block_size, n);

  /*--------------------Communication types for MPI--------------------------*/
  MPI_Datatype Block_a;
  MPI_Datatype Block_a_local;
  MPI_Datatype Block_b;
  MPI_Type_vector(m , Block_size, k_a, MPI_DOUBLE, &Block_a);
  MPI_Type_vector(m , Block_size, Block_size, MPI_DOUBLE, &Block_a_local);
  MPI_Type_vector(Block_size, n, n, MPI_DOUBLE, &Block_b);
  MPI_Type_commit(&Block_a);
  MPI_Type_commit(&Block_a_local);
  MPI_Type_commit(&Block_b);
  /*-------------Communication types for MPI are configured------------------*/


  MPI_Barrier(my_world);
  get_time(&start_time);

  if( NB_groups > 1 ) {
    hpnla_bcast(a, m*k_a, MPI_DOUBLE, 0, group_line, bcast_algorithm);
  /*  err = MPI_Bcast(a, m*k_a, MPI_DOUBLE, 0, group_line);    
    if (err != MPI_SUCCESS) {
      perror("Error Bcast A\n");
      return -1;
    }
   */
    hpnla_bcast(b, n*k_a, MPI_DOUBLE, 0, group_line, bcast_algorithm);
   /* err = MPI_Bcast(b, n*k_b, MPI_DOUBLE, 0, group_line);
    if (err != MPI_SUCCESS) {
      perror("Error Bcast B\n");
      return -1;
    }
   */
    MPI_Barrier(my_world);
  }

  get_time(&end_time_intern);
  communication_time += get_timediff(&start_time,&end_time_intern);

  info_print(0, row, col, "group %zu NB_block: %zu, NB_groups %zu\n"
              ,group,NB_Block, NB_groups);
  info_print(0, row, col,
              "m %zu,  k_a %zu,  k_b %zu,  n %zu,  Block_size %zu, "
              "group*NB_Block/NB_groups %zu, "
              "(group+1)*NB_Block/NB_groups %zu, row %zu,  col %zu,"
              "size_row %zu,  size_col %zu\n" ,
              m, k_a, k_b, n, Block_size,
              group*NB_Block/NB_groups, (group+1)*NB_Block/NB_groups,
              row, col, size_row, size_col);


  Summa( gemm, a, b, c, k_a, n, n, m, k_a, k_b, n, Block_size,
          group*NB_Block/NB_groups, (group+1)*NB_Block/NB_groups,
          row, col, size_row, size_col, a_Summa, b_Summa, Block_a,
          Block_a_local, Block_b, row_comm, col_comm, 0,
          &summa_computation_t,
          &summa_communication_t,
          bcast_algorithm,
          distribution);

  /*-------------------------End Summa algorihtm-----------------------------*/

  MPI_Comm_rank(group_line, &myrank);

  MPI_Barrier(my_world);
  get_time(&start_time_reduce);

  if( NB_groups > 1 ) {
    // a gather is better?
    err = MPI_Reduce(c, res, m*n, MPI_DOUBLE, MPI_SUM, 0, group_line);

    if (err != MPI_SUCCESS) {
      perror("Error reduce A\n");
      return -1;
    }

  }else{
    double *swap= c;
    c = res;
    res=swap;
  }

  MPI_Barrier(my_world);
  get_time(&end_time_reduce);

  MPI_Barrier(my_world);
  get_time(&end_time);

  time = get_timediff(&start_time,&end_time);
  double reduce_time = get_timediff(&start_time_reduce,&end_time_reduce);

  fprintf(stdout, "size_row: %zu, size_col: %zu, row: %zu, col: %zu, NB_proc: %d,NB_groups: %zu, n_global: %d, n: %zu, Block_size: %zu, Communication: %le, reduce: %le, "
         "Stime: %le, Scomputation: %le, Scommunication: %le, bcast_algorithm: %d, is_cyclic: %d\n",
         size_row,
         size_col,
         row,
         col,
         NB_proc,
         NB_groups,
         n_global,
         n,
         Block_size,
         communication_time,
         reduce_time,
         time,
         summa_computation_t,
         summa_communication_t,
         bcast_algorithm,
         distribution);


  MPI_Barrier(my_world);

#if CHECK_25D
  if(myrank == 0)
    check_result(res, a, b, m, n, k_a, k_b, row, col, size_row, size_col);
#endif

  // close properly the pragram

  MPI_Type_free(&Block_a);
  MPI_Type_free(&Block_a_local);
  MPI_Type_free(&Block_b);

  SMPI_SHARED_FREE(a_Summa);
  SMPI_SHARED_FREE(b_Summa);


  SMPI_SHARED_FREE(a);
  SMPI_SHARED_FREE(b);

  if( NB_groups > 1 ) {
    SMPI_SHARED_FREE(c);
  }

  SMPI_SHARED_FREE(res);

  MPI_Barrier(world);
  MPI_Comm_free(&my_world);
  MPI_Comm_free(&group_comm);
  MPI_Comm_free(&group_line);
  MPI_Comm_free(&row_comm);
  MPI_Comm_free(&col_comm);

  return 0;
}



