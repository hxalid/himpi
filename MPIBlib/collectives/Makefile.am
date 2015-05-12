AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)

lib_LTLIBRARIES = libmpib_coll.la
libmpib_coll_la_SOURCES = \
	mpib_coll.h \
	mpib_basic_collectives.h \
	mpib_tree_collectives.h \
	libmpib_coll.h \
	libmpib_coll.c \
	sgv_flat.h \
	sg_flat.c \
	sgv_flat.c \
	comm_tree.hpp \
	tree_algorithms.hpp \
	br_tree_algorithms.hpp \
	sg_tree_algorithms.hpp \
	sgv_tree_algorithms.hpp \
	brsg_tree_builders.hpp \
	sgv_tree_builders.hpp \
	mpib_tree_collectives.cpp

libmpib_coll_la_LIBADD = \
	$(top_builddir)/benchmarks/libmpib.la

EXTRA_DIST = \
	tree_collectives.dot
