AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) -I$(top_builddir)/utils -I$(top_builder)/hmpi

bin_PROGRAMS = test_hbcast
test_hbcast_SOURCES = \
    test_hbcast.c \
    mpi_bcast_algs.c \
    mpi_bcast_algs.h
nodist_EXTRA_test_hbcast_SOURCES = dummy.cpp

python_PYTHON = hmpi_runner.py
    
test_hbcast_LDADD = \
  $(top_builddir)/tools/libhmpi_tools.la \
  $(top_builddir)/hmpi/libhmpi.la \
  $(top_builddir)/MPIBlib/benchmarks/libmpib.la