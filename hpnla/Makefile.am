AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) \
              -I$(top_builddir)/cblas_wrappers -I$(top_srcdir)/cblas_wrappers

bin_PROGRAMS = cblas_test gemm_block gemm_block_mpi MM_mpi MM_recursive Cannon LU_runner

#Cannon example
Cannon_SOURCES = \
  MM_cannon.c \
  Cannon.c \
  Cannon.h \
  Matrix_init.c \
  Matrix_init.h
nodist_EXTRA_Cannon_SOURCES = dummy.cpp
Cannon_LDADD = \
  $(top_builddir)/cblas_wrappers/libcblas_wrappers.la \
  $(top_builddir)/tools/libhmpi_tools.la


cblas_test_SOURCES = cblas_test.c
nodist_EXTRA_cblas_test_SOURCES = dummy.cpp
cblas_test_LDADD = \
  $(top_builddir)/cblas_wrappers/libcblas_wrappers.la \
  $(top_builddir)/tools/libhmpi_tools.la
 

#Block Matrix Multiplication example
gemm_block_SOURCES = \
                     gemm_block.c \
                     Matrix_init.h \
                     Matrix_init.c
nodist_EXTRA_gemm_block_SOURCES = dummy.cpp
gemm_block_LDADD = \
  $(top_builddir)/cblas_wrappers/libcblas_wrappers.la \
  $(top_builddir)/tools/libhmpi_tools.la

#Summa example
gemm_block_mpi_SOURCES = \
  gemm_block_mpi.c \
  Summa.c \
  Summa.h \
  Matrix_init.c \
  Matrix_init.h
  
nodist_EXTRA_gemm_block_mpi_SOURCES = dummy.cpp
gemm_block_mpi_LDADD = \
  $(top_builddir)/cblas_wrappers/libcblas_wrappers.la \
  $(top_builddir)/tools/libhmpi_tools.la \
  $(top_builddir)/communication/libhpnla_bcast.la

#Summa 2.5D example
MM_mpi_SOURCES = \
  MM_mpi.c \
  Summa.c \
  Summa.h \
  Matrix_init.c \
  Matrix_init.h \
  2.5D_MM.c \
  2.5D_MM.h
  
nodist_EXTRA_MM_mpi_SOURCES = dummy.cpp
MM_mpi_LDADD = \
  $(top_builddir)/cblas_wrappers/libcblas_wrappers.la \
  $(top_builddir)/tools/libhmpi_tools.la \
  $(top_builddir)/communication/libhpnla_bcast.la

#Hierarchical Matrix Multiplication example
MM_recursive_SOURCES = \
  MM_recursive.c \
  Summa.c \
  Summa.h \
  Matrix_init.c \
  Matrix_init.h \
  Recursive.c \
  Recursive.h
  
nodist_EXTRA_MM_recursive_SOURCES = dummy.cpp
MM_recursive_LDADD = \
  $(top_builddir)/cblas_wrappers/libcblas_wrappers.la \
  $(top_builddir)/tools/libhmpi_tools.la \
  $(top_builddir)/communication/libhpnla_bcast.la



# 2D LU factorization
LU_runner_SOURCES = \
  LU_runner.c \
  Matrix_init.c \
  Matrix_init.h \
  lu_factorization.c \
  lu_factorization.h

nodist_EXTRA_LU_runner_SOURCES = dummy.cpp
LU_runner_LDADD = \
  $(top_builddir)/cblas_wrappers/libcblas_wrappers.la \
  $(top_builddir)/tools/libhmpi_tools.la \
  $(top_builddir)/communication/libhpnla_bcast.la

