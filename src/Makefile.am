AM_CFLAGS = -std=c99 -D CM_MODPOLDIR=\"${pkgdatadir}\" -I$(top_srcdir)/lib

bin_PROGRAMS = cm classpol ecpp
if MPI
bin_PROGRAMS += ecpp-mpi
endif

cm_LDADD = $(top_builddir)/lib/libcm.la

cm_SOURCES = params.h params.c cm.c

classpol_LDADD = $(top_builddir)/lib/libcm.la

classpol_SOURCES = params.h params.c classpol.c

ecpp_LDADD = $(top_builddir)/lib/libcm.la

ecpp_SOURCES = params.h params.c ecpp.c

if MPI
ecpp_mpi_LDADD = $(top_builddir)/lib/libcm_mpi.la

ecpp_mpi_SOURCES = params.h params.c ecpp-mpi.c

ecpp_mpi_CFLAGS = $(AM_CFLAGS) $(MPI_CFLAGS) -D WITH_MPI

ecpp_mpi_LDFLAGS = $(MPI_CLDFLAGS)
endif

