AM_CFLAGS=-std=c99

include_HEADERS = cm_common.h cm_arith.h cm_class.h

lib_LTLIBRARIES = libcm.la

libcm_la_LDFLAGS = -version-info 0:0:0

libcm_la_SOURCES = cm_common.h cm_arith.h cm_common-impl.h \
   timer.c file.c nt.c qdev.c modular.c modpol.c fem.c \
   mpzx.c classgroup.c modclass.c curve.c class.c pari.c
