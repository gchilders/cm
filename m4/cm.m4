#
# SYNOPSIS
#
#
# MPC_C_CHECK_FLAG([FLAG,ACCUMULATOR])
#
# DESCRIPTION
#
# Checks if the C compiler accepts the flag FLAG
# If yes, adds it to CFLAGS.

AC_DEFUN([MPC_C_CHECK_FLAG], [
   AX_C_CHECK_FLAG($1,,,[CFLAGS="$CFLAGS $1"])
])


#
# SYNOPSIS
#
#
# MPC_C_CHECK_WARNINGCFLAGS
#
# DESCRIPTION
#
# For development version only: Checks if gcc accepts warning flags.
# Adds accepted ones to CFLAGS.
#
AC_DEFUN([MPC_C_CHECK_WARNINGCFLAGS], [
  AC_REQUIRE([AC_PROG_GREP])
  if echo $VERSION | grep -c dev >/dev/null 2>&1 ; then
    if test "x$GCC" = "xyes" -a "x$compiler" != "xicc" -a "x$compiler" != "xg++"; then
      case $host in
         *darwin*) ;;
         *) MPC_C_CHECK_FLAG(-D_FORTIFY_SOURCE=2,$1) ;;
      esac
      MPC_C_CHECK_FLAG(-g)
      MPC_C_CHECK_FLAG(-std=c99)
      MPC_C_CHECK_FLAG(-pedantic)
      MPC_C_CHECK_FLAG(-Wno-long-long)
      MPC_C_CHECK_FLAG(-Wall)
      MPC_C_CHECK_FLAG(-Wextra)
      MPC_C_CHECK_FLAG(-Werror)
      MPC_C_CHECK_FLAG(-Wdeclaration-after-statement)
      MPC_C_CHECK_FLAG(-Wshadow)
      MPC_C_CHECK_FLAG(-Wstrict-prototypes)
      MPC_C_CHECK_FLAG(-Wmissing-prototypes)
      MPC_C_CHECK_FLAG(-Wno-unused-value)
    fi
  fi
])
