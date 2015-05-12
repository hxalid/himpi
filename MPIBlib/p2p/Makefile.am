AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) 

lib_LTLIBRARIES = libmpib_p2p.la
libmpib_p2p_la_SOURCES = \
	libmpib_p2p.h \
	libmpib_p2p.c \
	mpib_p2p.h \
	mpib_p2p.c \
	mpib_p2p_sg.c

libmpib_p2p_la_LIBADD = $(top_builddir)/benchmarks/libmpib.la
