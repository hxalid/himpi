AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)

include_HEADERS = \
        hpnla_debug.h \
        hpnla_conf.h \
        hpnla_bind.h \
        hpnla_timer.h \
        utils.h \
        optimal_groups.h

noinst_LTLIBRARIES = libhmpi_tools.la
libhmpi_tools_la_SOURCES = \
        hpnla_conf.c \
        hpnla_bind.c \
        hpnla_timer.c \
        utils.c \
        optimal_groups.c
