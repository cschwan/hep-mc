# initialize configure
AC_INIT([hep-mc], [0.2.1.9999], [https://github.com/cschwan/hep-mc/issues],
	[hep-mc], [https://github.com/cschwan/hep-mc])

# initialize automake
AM_INIT_AUTOMAKE([foreign -Wall dist-bzip2 no-dist-gzip nostdinc])

# make automake silent if possible
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

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

# enable building of examples
AC_ARG_ENABLE([examples], AS_HELP_STRING([--enable-examples],
	[Enable building of example programs]))

AM_CONDITIONAL([HAVE_EXAMPLES], [test "x$enable_examples" = "xyes"])

# some examples need MPI
AS_IF([test "x$enable_examples" = "xyes"], [LX_FIND_MPI])

# add possibility to generate API documentation with Doxygen
AC_ARG_ENABLE([doxygen], AS_HELP_STRING([--enable-doxygen],
	[Enable generation of Doxygen documentation]))

# variable pointing to the doxygen binary
AC_ARG_VAR([DOXYGEN], [Doxygen documentation generator])

# check for the doxygen binary
AS_IF([test "x$enable_doxygen" = "xyes"], [
	AC_CHECK_PROGS([DOXYGEN], [doxygen])
	AS_IF([test -z "$DOXYGEN"],
		AC_MSG_ERROR([Doxygen not found but requested]))
])

# if doxygen is requested, define HAVE_DOXYGEN and configure Doxyfile.in
AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])

# check if the user wants to test this library
AC_ARG_ENABLE([tests], AS_HELP_STRING([--enable-tests],
	[Enable unit tests (requires Boost Unit Test Framework)]))

# check for header and library
AS_IF([test "x$enable_tests" = "xyes"], [
	# save content of LIBS variable
	previous_libs="${LIBS}"

	AC_CHECK_HEADERS([boost/mpl/list.hpp], [],
		AC_MSG_ERROR([header not found/usable]))
	AC_CHECK_HEADERS([boost/test/unit_test.hpp], [],
		AC_MSG_ERROR([header not found/usable]))
	AC_SEARCH_LIBS([_ZN5boost9unit_test14unit_test_mainEPFbvEiPPc],
		[boost_unit_test_framework], [], AC_MSG_ERROR([library not found]))

	# restore LIBS variable - boost libs are only needed in unit tests
	LIBS="${previous_libs}"
])

# if tests are requested, define HAVE_TESTS
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
