AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) -I$(top_builddir)/hmpi

include_HEADERS = hmpi.h

noinst_LTLIBRARIES = libhmpi.la
libhmpi_la_SOURCES = \
    	mpi_bcast_algs.c \
    	mpi_bcast_algs.h \
        hmpi.h \
		hmpi_init.c \
		hmpi_finalize.c \
		hallreduce.c \
        hbcast.c \
        hreduce.c \
        hgather.c \
        hscatter.c 
        
libhmpi_la_LIBADD = \
  $(top_builddir)/tools/libhmpi_tools.la