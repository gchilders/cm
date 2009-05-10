AC_PREREQ(2.61)
AC_INIT(cm, 0.1-dev, enge@lix.polytechnique.fr)
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([lib/common/cm_common.h])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([1.9])
AM_MAINTAINER_MODE

# Extra arguments to configure
AC_ARG_WITH([gmp],
            [AC_HELP_STRING([--with-gmp=DIR],
                           [GMP install directory])],
            [
            CPPFLAGS="$CPPFLAGS -I$withval/include"
            LDFLAGS="$LDFLAGS -L$withval/lib"
            ])
AC_ARG_WITH([mpfr],
            [AC_HELP_STRING([--with-mpfr=DIR],
                           [MPFR install directory])],
            [
            CPPFLAGS="$CPPFLAGS -I$withval/include"
            LDFLAGS="$LDFLAGS -L$withval/lib"
            ])
AC_ARG_WITH([mpc],
            [AC_HELP_STRING([--with-mpc=DIR],
                           [MPC install directory])],
            [
            CPPFLAGS="$CPPFLAGS -I$withval/include"
            LDFLAGS="$LDFLAGS -L$withval/lib"
            ])
AC_ARG_WITH([mpfrcx],
            [AC_HELP_STRING([--with-mpfrcx=DIR],
                           [MPFRCX install directory])],
            [
            CPPFLAGS="$CPPFLAGS -I$withval/include"
            LDFLAGS="$LDFLAGS -L$withval/lib"
            ])
AC_ARG_WITH([ntl],
            [AC_HELP_STRING([--with-ntl=DIR],
                           [NTL install directory])],
            [
            CPPFLAGS="$CPPFLAGS -I$withval/include"
            LDFLAGS="$LDFLAGS -L$withval/lib"
            ])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([inttypes.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_FUNCS([sqrt])

# Checks if compiler accepts warning flags.
CM_PROG_CC_WARNINGCFLAGS([WARNINGCFLAGS])

# Checks for libraries.
AC_CHECK_LIB([z], [gzopen])
AC_CHECK_LIB([gmp], [__gmpz_init],
             [LIBS="-lgmp $LIBS"],
             [AC_MSG_ERROR([libgmp not found or uses a different ABI.])])
AC_MSG_CHECKING(for MPFR)
LIBS="-lmpfr $LIBS"
AC_LINK_IFELSE(
        [AC_LANG_PROGRAM(
                [[#include "mpfr.h"]],
                [[mpfr_t x;  mpfr_init(x) ; mpfr_clear(x);]]
        )],
        [AC_MSG_RESULT(yes)],
        [
        AC_MSG_RESULT(no)
        AC_MSG_ERROR([libmpfr not found or uses a different ABI.])
        ])
AC_MSG_CHECKING(for MPC)
LIBS="-lmpc $LIBS"
AC_LINK_IFELSE(
        [AC_LANG_PROGRAM(
                [[#include "mpc.h"]],
                [[mpc_t x;  mpc_init2(x, 2) ; mpc_clear(x);]]
        )],
        [AC_MSG_RESULT(yes)],
        [
        AC_MSG_RESULT(no)
        AC_MSG_ERROR([libmpc not found or uses a different ABI.])
        ])
AC_MSG_CHECKING(for MPFRCX)
LIBS="-lmpfrcx $LIBS"
AC_LINK_IFELSE(
        [AC_LANG_PROGRAM(
                [[#include "mpfrcx.h"]],
                [[mpcx_t x;  mpcx_init(x, 2, 2) ; mpcx_clear(x);]]
        )],
        [AC_MSG_RESULT(yes)],
        [
        AC_MSG_RESULT(no)
        AC_MSG_ERROR([libmpfrcx not found or uses a different ABI.])
        ])
AC_LANG_PUSH([C++])
AC_MSG_CHECKING(for NTL)
LIBS="-lntl $LIBS"
AC_LINK_IFELSE(
        [AC_LANG_PROGRAM(
                [[#include "NTL/ZZ.h"]],
                [[NTL::ZZ x;]]
        )],
        [AC_MSG_RESULT(yes)],
        [
        AC_MSG_RESULT(no)
        AC_MSG_ERROR([libntl not found or uses a different ABI.])
        ])
AC_LANG_POP([C++])

AC_CONFIG_FILES([Makefile lib/Makefile lib/mpfp/Makefile lib/mpfpx/Makefile lib/common/Makefile lib/class/Makefile data/Makefile data/df/Makefile src/Makefile test/Makefile])
AC_OUTPUT
