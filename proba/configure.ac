AC_INIT([frmw_proba], [0.1], [mjoaojr@gmail.com])

# Detect the canonical host and target build environment
# see e.g. http://cvs.haskell.org/darcs/ghc/configure.ac for howto use
# AC_CANONICAL_HOST or AC_CANONICAL_TARGET?
# AC_CANONICAL_SYSTEM will also run AC_CANONICAL_HOST and allows to use --target
AC_CANONICAL_SYSTEM
# DEFINES here are going to config.h and have to be the same as in clustalx.pro
#echo "host_os $host_os"
case $target_os in
   *linux*|*solaris*|*bsd*)
        AC_DEFINE([OS_UNIX], [1], [Host OS])
        ;;
   *darwin*)
        AC_DEFINE([OS_MAC], [1], [Host OS])
        AC_DEFINE([OS_UNIX], [1], [Host OS])
        ;;
   *cygwin*|*mingw*)
        AC_DEFINE([OS_WINDOWS], [1], [Host OS])
        ;;
   *)
        AC_MSG_WARN(["Couldn't determine target operating system. Using Unix"]);
        AC_DEFINE([OS_UNIX], [1], [Host OS])
        ;;
esac
echo "building for $target_os"

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_PROG_CXX
AC_PROG_CC
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT
