AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) -I$(top_builddir)/hmpi

include_HEADERS = hmpi.h

noinst_LTLIBRARIES = libhmpi.la
libhmpi_la_SOURCES = \
        hmpi.h \
		hmpi_init.c \
		hmpi_finalize.c \
		hallreduce.c \
        hbcast.c \
        hgather.c \
        hreduce.c \
        hscatter.c