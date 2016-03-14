AM_CFLAGS=-std=c99 -I$(top_srcdir)/lib/mpfp -I$(top_srcdir)/lib/mpfpx -I$(top_srcdir)/lib/common
AM_CXXFLAGS=@WARNINGCFLAGS@ -I$(top_srcdir)/lib/mpfp -I$(top_srcdir)/lib/mpfpx -I$(top_srcdir)/lib/common

include_HEADERS = cm_class.h

lib_LTLIBRARIES = libcm_class.la

libcm_class_la_LDFLAGS = -version-info 1:0:0

libcm_class_la_SOURCES = cm_class.h cm_class-impl.h \
   classgroup.c modclass.c curve.c class.c pari.c
