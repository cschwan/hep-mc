# AC_CONFIG_MACRO_DIR needs >=2.58
AC_PREREQ([2.58])

# initialize configure
AC_INIT([hep-mc], [0.2.2.9999], [https://github.com/cschwan/hep-mc/issues],
	[hep-mc], [https://github.com/cschwan/hep-mc])

# put autotools files into a separate directory
AC_CONFIG_AUX_DIR([build-aux])

# initialize automake
AM_INIT_AUTOMAKE([foreign -Wall dist-bzip2 no-dist-gzip nostdinc])

# put aclocal macros into m4
AC_CONFIG_MACRO_DIR([m4])

# C++ is our language
AC_LANG([C++])

# check for C++11 conforming compiler (autoconf-archive macro)
AX_CXX_COMPILE_STDCXX_11([noext])

# point to a source file in order to let ./configure find itself
AC_CONFIG_SRCDIR([include/hep/mc.hpp])

# do not pollute command line with compiler definitions but instead write them
# into the specified header
AC_CONFIG_HEADERS([config.hpp])

# MPI
AX_MPI_FLAGS_INIT
AC_ARG_ENABLE([mpi], AS_HELP_STRING([--enable-mpi], [Enable use of the Message Passing Interface]))
AM_CONDITIONAL([HAVE_MPI], [test "x$enable_mpi" = "xyes"])
AS_IF([test "x$enable_mpi" = "xyes"], [AX_MPI_FLAGS])

# EXAMPLES
AC_ARG_ENABLE([examples], AS_HELP_STRING([--enable-examples],
	[Enable building of example programs]))
AM_CONDITIONAL([HAVE_EXAMPLES], [test "x$enable_examples" = "xyes"])

# DOXYGEN
AC_ARG_ENABLE([doxygen], AS_HELP_STRING([--enable-doxygen],
	[Enable generation of Doxygen documentation]))
AC_ARG_VAR([DOXYGEN], [Doxygen documentation generator])
AS_IF([test "x$enable_doxygen" = "xyes"], [
	AC_CHECK_PROGS([DOXYGEN], [doxygen])
	AS_IF([test -z "$DOXYGEN"], AC_MSG_ERROR([Doxygen not found but requested]))
])
AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])

# UNIT TESTS
AC_ARG_ENABLE([tests], AS_HELP_STRING([--enable-tests], [Enable unit tests]))
AS_IF([test "x$enable_tests" = "xyes"], [
	# fix -Wextra-portability warning
	m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

	# initialize libtool
	LT_INIT
])
AM_CONDITIONAL([HAVE_TESTS], [test "x$enable_tests" = "xyes"])

# configure Makefiles and Doxyfile
AC_CONFIG_FILES([
	Makefile
	doc/Doxyfile
	doc/Makefile
	examples/Makefile
	include/Makefile
	tests/Makefile
])

# write configure script
AC_OUTPUT
