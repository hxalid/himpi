AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)

noinst_PROGRAMS = substitute sgv tree_builders sg-based-p2p p2p-eager fragment-omp-for

substitute_SOURCES = substitute.c
substitute_CPPFLAGS = -I$(top_srcdir)/collectives
substitute_LDADD = $(top_builddir)/collectives/libmpib_coll.la

sgv_SOURCES = sgv.c
sgv_LDADD = $(top_builddir)/collectives/libmpib_coll.la

tree_builders_SOURCES = tree_builders.cpp

sg_based_p2p_SOURCES = sg-based-p2p.c
sg_based_p2p_LDADD = $(top_builddir)/p2p/libmpib_p2p.la

p2p_eager_SOURCES = p2p-eager.c

fragment_omp_for_SOURCES = fragment-omp-for.c
fragment_omp_for_CFLAGS = -fopenmp 
