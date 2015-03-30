AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) -I$(top_builddir)/utils -I$(top_builder)/hmpi

bin_PROGRAMS = test_hbcast
test_hbcast_SOURCES = \
    test_hbcast.c \
    mpi_bcast_algs.c \
    mpi_bcast_algs.h 

    
test_hbcast_LDADD = \
  $(top_builddir)/tools/libhpnla_tools.la \
  $(top_builddir)/hmpi/libhmpi.la