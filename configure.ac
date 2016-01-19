# AC_CONFIG_MACRO_DIR needs >=2.58
AC_PREREQ([2.58])

# initialize configure
AC_INIT([hep-mc], [0.5], [https://github.com/cschwan/hep-mc/issues],
	[hep-mc], [https://github.com/cschwan/hep-mc])

# put autotools files into a separate directory
AC_CONFIG_AUX_DIR([build-aux])

# initialize automake
AM_INIT_AUTOMAKE([1.11 foreign -Wall dist-bzip2 no-dist-gzip nostdinc])

# put aclocal macros into m4
AC_CONFIG_MACRO_DIR([m4])

# C++ is our language
AC_LANG([C++])

# unit tests need convenience libraries
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_RANLIB

# check for C++11 conforming compiler (autoconf-archive macro)
AX_CXX_COMPILE_STDCXX([11], [], [mandatory])

# point to a source file in order to let ./configure find itself
AC_CONFIG_SRCDIR([include/hep/mc.hpp])

# don't pollute the command line with compiler defs, write them into a header that we dont need
AC_CONFIG_HEADERS([config.hpp])

# MPI
AC_ARG_VAR(MPICXX, [MPI C++ compiler command])
AC_ARG_ENABLE([mpi], AS_HELP_STRING([--enable-mpi], [Enable use of the Message Passing Interface]))
AM_CONDITIONAL([HAVE_MPI], [test "x$enable_mpi" = "xyes"])
AS_IF([test "x$enable_mpi" = "xyes"], [
	AC_PATH_PROGS([MPICXX], [mpicxx mpiCC mpic++ mpig++ mpiicpc mpipgCC mpixlC])
	AS_IF([test x"$MPICXX" == x""], AC_MSG_ERROR([MPI requested but no appropriate compiler found.]))
	AX_MPI_FLAGS([MPICXX])
])

# EXAMPLES
AC_ARG_ENABLE([examples], AS_HELP_STRING([--enable-examples],
	[Enable building of example programs]))
AM_CONDITIONAL([HAVE_EXAMPLES], [test "x$enable_examples" = "xyes"])

# DOXYGEN
AC_ARG_ENABLE([doxygen], AS_HELP_STRING([--enable-doxygen],
	[Enable generation of Doxygen documentation]))
AC_ARG_VAR([DOXYGEN], [Doxygen documentation generator])
AS_IF([test "x$enable_doxygen" = "xyes"], [
	AC_PATH_PROGS([DOXYGEN], [doxygen])
	AS_IF([test -z "$DOXYGEN"], AC_MSG_ERROR([Doxygen not found but requested]))
])
AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])

# configure Makefiles
AC_CONFIG_FILES([
	Makefile
	doc/Makefile
	examples/Makefile
	include/Makefile
	tests/Makefile
])

# configure remaining files
AC_CONFIG_FILES([doc/Doxyfile])
AC_CONFIG_FILES([tests/test_mpi_numerical_results.sh],
	[chmod +x tests/test_mpi_numerical_results.sh])

# write configure script
AC_OUTPUT
