AM_CFLAGS = -std=c99 -D CM_MODPOLDIR=\"${pkgdatadir}\" -I$(top_srcdir)/lib/mpfpx -I$(top_srcdir)/lib/common -I$(top_srcdir)/lib/class

bin_PROGRAMS = cm classpol

cm_LDADD = $(top_builddir)/lib/class/libcm_class.la $(top_builddir)/lib/common/libcm_common.la $(top_builddir)/lib/mpfpx/libmpfpx.la

cm_SOURCES = params.h params.c cm.c

classpol_LDADD = $(top_builddir)/lib/class/libcm_class.la $(top_builddir)/lib/common/libcm_common.la $(top_builddir)/lib/mpfpx/libmpfpx.la

classpol_SOURCES = params.h params.c classpol.c
