AM_CFLAGS = -std=c99 -D CM_MODPOLDIR=\"${datadir}\" @WARNINGCFLAGS@ -I$(top_srcdir)/lib/mpfp -I$(top_srcdir)/lib/mpfpx -I$(top_srcdir)/lib/common -I$(top_srcdir)/lib/class
AM_CXXFLAGS = -D CM_MODPOLDIR=\"${datadir}\" @WARNINGCFLAGS@ -I$(top_srcdir)/lib/mpfp -I$(top_srcdir)/lib/mpfpx -I$(top_srcdir)/lib/common -I$(top_srcdir)/lib/class
AM_LDFLAGS = -all-static

bin_PROGRAMS = cm classpol

cm_LDADD = $(top_builddir)/lib/class/libcm_class.la $(top_builddir)/lib/common/libcm_common.la $(top_builddir)/lib/mpfpx/libmpfpx.la $(top_builddir)/lib/mpfp/libmpfp.la

cm_SOURCES = params.h params.c cm.c

classpol_LDADD = $(top_builddir)/lib/class/libcm_class.la $(top_builddir)/lib/common/libcm_common.la $(top_builddir)/lib/mpfpx/libmpfpx.la $(top_builddir)/lib/mpfp/libmpfp.la

classpol_SOURCES = params.h params.c classpol.c
