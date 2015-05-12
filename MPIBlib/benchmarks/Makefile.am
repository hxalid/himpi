AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)

noinst_LTLIBRARIES = libmpib.la
libmpib_la_SOURCES = \
	mpib.h \
	mpib_measurement.h \
	mpib_measurement.c \
	mpib_defs.h \
	mpib_benchmarks.h \
	mpib_measure_bcast.c \
	mpib_p2p_benchmarks.h \
	mpib_p2p_benchmarks.c \
	mpib_coll_benchmarks.h \
	mpib_measure_max.c \
	mpib_measure_root.c \
	mpib_measure_global.c \
	mpib_coll_containers.h \
	mpib_coll_containers.hpp \
	mpib_coll_containers.cpp \
	mpib_getopt.h \
	mpib_getopt.c \
	mpib_shared.h \
	mpib_output.h \
	mpib_output.c \
	mpib_utilities.c \
	mpib_utilities.h

EXTRA_DIST = \
	benchmark_msgset.dot \
	benchmark.dot \
	coll_containers.dot \
	modules.dot
