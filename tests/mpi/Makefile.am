AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) -I$(top_builddir)/utils -I$(top_builder)/himpi

bin_PROGRAMS = test_hbcast
test_hbcast_SOURCES = \
    test_hbcast.c 
nodist_EXTRA_test_hbcast_SOURCES = dummy.cpp

python_PYTHON = hmpi_runner.py
    
test_hbcast_LDADD = \
  $(top_builddir)/tools/libhimpi_tools.la \
  $(top_builddir)/himpi/libhimpi.la \
  $(top_builddir)/MPIBlib/benchmarks/libmpib.la