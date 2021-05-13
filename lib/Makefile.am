AM_CFLAGS=-std=c99

include_HEADERS = cm.h cm-arith.h

lib_LTLIBRARIES = libcm.la

libcm_la_LDFLAGS = -version-info 0:0:0

libcm_la_SOURCES = cm.h cm-arith.h cm-impl.h \
   timer.c file.c nt.c qdev.c modular.c modpol.c mpzx.c classgroup.c \
   modclass.c curve.c param.c class.c jmodp.c pari.c
