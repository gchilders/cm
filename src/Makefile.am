AM_CFLAGS = -std=c99 -D CM_MODPOLDIR=\"${pkgdatadir}\" -I$(top_srcdir)/lib

bin_PROGRAMS = cm classpol

cm_LDADD = $(top_builddir)/lib/libcm.la

cm_SOURCES = params.h params.c cm.c

classpol_LDADD = $(top_builddir)/lib/libcm.la

classpol_SOURCES = params.h params.c classpol.c
