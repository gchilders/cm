AM_CFLAGS=-std=c99 -D CM_MODPOLDIR=\"${top_srcdir}/data\" -I$(top_srcdir)/lib/mpfpx -I$(top_srcdir)/lib/common -I$(top_srcdir)/lib/class

AM_LDFLAGS = -no-install

LDADD = $(top_builddir)/lib/class/libcm_class.la $(top_builddir)/lib/common/libcm_common.la $(top_builddir)/lib/mpfpx/libmpfpx.la

check_PROGRAMS = tcm

TESTS_ENVIRONMENT = $(VALGRIND)
TESTS = $(check_PROGRAMS)
