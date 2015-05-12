AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)

bin_PROGRAMS = \
	p2p \
	collective \
	collective_test \
	generate_factors \
	p2p_forward \
	p2p_pairs

p2p_SOURCES = p2p.c
nodist_EXTRA_p2p_SOURCES = dummy.cpp
p2p_LDADD = $(top_builddir)/p2p/libmpib_p2p.la

collective_SOURCES = collective.c
nodist_EXTRA_collective_SOURCES = dummy.cpp
collective_LDADD = $(top_builddir)/collectives/libmpib_coll.la

collective_test_SOURCES = \
	collective_test.c \
	mpib_test.c \
	mpib_test.h
collective_test_LDADD = $(top_builddir)/benchmarks/libmpib.la

generate_factors_SOURCES = generate_factors.c
nodist_EXTRA_collective_test_SOURCES = dummy.cpp

p2p_forward_SOURCES = p2p_forward.c
p2p_forward_LDADD = $(top_builddir)/p2p/libmpib_p2p.la

p2p_pairs_SOURCES = p2p_pairs.c
p2p_pairs_LDADD = $(top_builddir)/benchmarks/libmpib.la

EXTRA_DIST = \
	p2p.plot \
	collective.plot
