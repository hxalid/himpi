AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) -I$(top_builddir)/utils

bin_PROGRAMS = test_hbcast
test_hbcast_SOURCES = \
    test_hbcast.c \
    hbcast_old.c \
    hbcast_old.h \
    mpi_bcast_algs.c \
    mpi_bcast_algs.h 

    
test_hbcast_LDADD = \
  $(top_builddir)/tools/libhpnla_tools.la