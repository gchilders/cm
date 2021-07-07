/*
configure.ac

Copyright (C) 2009, 2010, 2011, 2012, 2013, 2015, 2016, 2018, 2021 Andreas Enge

This file is part of CM.

CM is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the license, or (at your
option) any later version.

CM is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with CM; see the file COPYING. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

AC_PREREQ([2.69])
AC_INIT([cm],[0.4.0dev],[andreas.enge@inria.fr])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([lib/cm.h])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

AC_CANONICAL_HOST
AC_CONFIG_MACRO_DIR([m4])

# Extra arguments to configure
AC_ARG_WITH([gmp],
            [AS_HELP_STRING([--with-gmp=DIR],
                           [GMP install directory])],
            [
            CPPFLAGS="$CPPFLAGS -I$withval/include"
            LDFLAGS="$LDFLAGS -L$withval/lib"
            ])
AC_ARG_WITH([mpfr],
            [AS_HELP_STRING([--with-mpfr=DIR],
                           [MPFR install directory])],
            [
            CPPFLAGS="$CPPFLAGS -I$withval/include"
            LDFLAGS="$LDFLAGS -L$withval/lib"
            ])
AC_ARG_WITH([mpc],
            [AS_HELP_STRING([--with-mpc=DIR],
                           [MPC install directory])],
            [
            CPPFLAGS="$CPPFLAGS -I$withval/include"
            LDFLAGS="$LDFLAGS -L$withval/lib"
            ])
AC_ARG_WITH([mpfrcx],
            [AS_HELP_STRING([--with-mpfrcx=DIR],
                           [MPFRCX install directory])],
            [
            CPPFLAGS="$CPPFLAGS -I$withval/include"
            LDFLAGS="$LDFLAGS -L$withval/lib"
            ])
AC_ARG_WITH([pari],
            [AS_HELP_STRING([--with-pari=DIR],
                           [PARI install directory])],
            [
            CPPFLAGS="$CPPFLAGS -I$withval/include"
            LDFLAGS="$LDFLAGS -L$withval/lib -Wl,-rpath -Wl,$withval/lib"
            ])
AC_ARG_ENABLE([valgrind-tests],
             [AS_HELP_STRING([--enable-valgrind-tests],
                             [run checks through valgrind (default = no)])],
             [case $enableval in
                 yes) gl_VALGRIND_TESTS ;;
                 no)  ;;
                 *)   AC_MSG_ERROR([Bad value for --enable-valgrind-tests: Use yes or no]) ;;
              esac
             ]
             )

# Checks for programs.
AC_PROG_CC
AC_LANG(C)

# Set up LibTool
LT_INIT

# Check GMP Header
AC_MSG_CHECKING(for gmp.h)
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include "gmp.h"
]])],[AC_MSG_RESULT(yes)],[
 AC_MSG_RESULT(no)
 AC_MSG_ERROR([gmp.h cannot be found or is unusable.])
])

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_LIB([m], [sqrt], , AC_MSG_ERROR([libm not found.]))
AC_CHECK_LIB([z], [gzopen], , AC_MSG_ERROR([libz not found.]))
AC_CHECK_FUNCS([sqrt])

# Checks if compiler accepts warning flags.
MPC_C_CHECK_WARNINGCFLAGS

# Gmp
AC_CHECK_LIB([gmp], [__gmpz_init],
             [LIBS="-lgmp $LIBS"],
             [AC_MSG_ERROR([libgmp not found or uses a different ABI.])])
AC_MSG_CHECKING(for recent GMP)
AC_COMPILE_IFELSE(
        [AC_LANG_SOURCE(
        [[
#include "gmp.h"
#if (__GNU_MP_VERSION*100 + __GNU_MP_VERSION_MINOR*10 + __GNU_MP_VERSION_PATCHLEVEL < 432)
# error "Minimal GMP version is 4.3.2"
error
#endif
        ]])],
        [AC_MSG_RESULT(yes)],
        [
        AC_MSG_RESULT(no)
        AC_MSG_ERROR([GMP version >= 4.3.2 required])
        ])

# Mpfr
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
AC_MSG_CHECKING(for recent MPFR)
AC_COMPILE_IFELSE(
        [AC_LANG_SOURCE(
        [[
#include "mpfr.h"
#if (MPFR_VERSION < MPFR_VERSION_NUM (3,0,0))
# error "Minimal MPFR version is 3.0.0"
error
#endif
        ]])],
        [AC_MSG_RESULT(yes)],
        [
         AC_MSG_RESULT(no)
         AC_MSG_ERROR([MPFR version >= 3.0.0 required])
         ])

# Mpc
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
AC_MSG_CHECKING(for recent MPC)
AC_COMPILE_IFELSE(
        [AC_LANG_SOURCE(
        [[
#include "mpc.h"
#if (MPC_VERSION < MPC_VERSION_NUM (1,0,0))
# error "Minimal MPC version is 1.0.0"
error
#endif
        ]])],
        [AC_MSG_RESULT(yes)],
        [
         AC_MSG_RESULT(no)
         AC_MSG_ERROR([MPC version >= 1.0.0 required])
         ])

# Mpfrcx
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
AC_MSG_CHECKING(for recent MPFRCX)
AC_COMPILE_IFELSE(
        [AC_LANG_SOURCE(
        [[
#include "mpfrcx.h"
#if (10000*MPFRCX_VERSION_MAJOR + 100*MPFRCX_VERSION_MINOR + MPFRCX_VERSION_PATCHLEVEL < 603)
# error "Minimal MPFRCX version is 0.6.3"
error
#endif
        ]])],
        [AC_MSG_RESULT(yes)],
        [
         AC_MSG_RESULT(no)
         AC_MSG_ERROR([MPFRCX version >= 0.6.3 required])
         ])

# Pari
AC_MSG_CHECKING(for PARI)
LIBS="-lpari $LIBS"
AC_LINK_IFELSE(
        [AC_LANG_PROGRAM(
                [[#include "pari/pari.h"]],
                [[pari_init (100, 0); pari_close ();]]
        )],
        [AC_MSG_RESULT(yes)],
        [
        AC_MSG_RESULT(no)
        AC_MSG_ERROR([libpari not found or uses a different ABI.])
        ])
AC_MSG_CHECKING(for recent PARI)
AC_COMPILE_IFELSE(
        [AC_LANG_SOURCE(
        [[
#include "pari/pari.h"
#if PARI_VERSION_CODE < PARI_VERSION (2,11,0)
# error "Minimal PARI version is 2.11.0"
error
#endif
        ]])],
        [AC_MSG_RESULT(yes)],
        [
         AC_MSG_RESULT(no)
         AC_MSG_ERROR([PARI version >= 2.11.0 required])
         ])

# Checks for gcc version, result in variables CM_GCC etc. inside config.h
AX_GCC_VERSION
AC_DEFINE_UNQUOTED([CM_GCC], ["$GCC"], [Gcc yes or no])
AC_DEFINE_UNQUOTED([CM_CC], ["$CC"], [C compiler])
AC_DEFINE_UNQUOTED([CM_GCC_VERSION], ["$GCC_VERSION"], [Version of gcc])


AC_CONFIG_FILES([Makefile doc/Makefile lib/Makefile data/Makefile data/af/Makefile data/df/Makefile data/mf/Makefile src/Makefile tests/Makefile])
AC_OUTPUT
