LOG_COMPILER = $(LOG_VALGRIND)

AM_CFLAGS=-std=c99 -D CM_MODPOLDIR=\"${top_srcdir}/data\" -I$(top_srcdir)/lib

AM_LDFLAGS = -no-install

LDADD = $(top_builddir)/lib/libcm.la

check_PROGRAMS = tcm ttwist tecpp

TESTS_ENVIRONMENT = $(VALGRIND)
TESTS = $(check_PROGRAMS)
