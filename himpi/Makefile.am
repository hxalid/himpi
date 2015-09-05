AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) -I$(top_builddir)/himpi

include_HEADERS = himpi.h

noinst_LTLIBRARIES = libhimpi.la
libhimpi_la_SOURCES = \
    	mpi_bcast_algs.c \
    	mpi_bcast_algs.h \
        himpi.h \
		himpi_init.c \
		himpi_finalize.c \
		hallreduce.c \
        hbcast.c \
        hreduce.c \
        hgather.c \
        hscatter.c 
        
libhimpi_la_LIBADD = \
  $(top_builddir)/tools/libhimpi_tools.la