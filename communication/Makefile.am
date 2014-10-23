AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) -I$(top_builddir)/communication


bin_PROGRAMS = hbcast_runner


include_HEADERS = \
	hpnla_bcast.h \
        hbcast.h

noinst_LTLIBRARIES = \
	libhpnla_bcast.la \
	libhbcast.la

libhpnla_bcast_la_SOURCES = \
	hpnla_bcast.c \
	hpnla_bcast.h 

libhbcast_la_SOURCES = \
        hbcast.c \
        hbcast.h


# hbast_runner
hbcast_runner_SOURCES = \
  hbcast_runner.c

nodist_EXTRA_hbcast_runner_SOURCES = dummy.cpp
hbcast_runner_LDADD = \
  $(top_builddir)/tools/libhpnla_tools.la \
  $(top_builddir)/communication/libhpnla_bcast.la \
  $(top_builddir)/communication/libhbcast.la

