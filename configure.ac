AC_PREREQ([2.63b])

# check if git is available and current directory contains repository
m4_esyscmd_s([git status])

# if a repository is available use git to create files that are needed before
# AC_INIT is called:
ifelse(sysval,
	# VERSION file, ...
	m4_esyscmd_s([echo `git describe 2>/dev/null` | tee VERSION])

	# AUTHORS file, ...
	m4_esyscmd_s([git log --format='%aN <%aE>' | uniq > AUTHORS])

	# ChangeLog file and ...
	m4_esyscmd_s([
		git log --date=short \
			--format='%ad  %aN  <%aE> %n    %w(80,0,4)%B' > ChangeLog
	])

	# NEWS file
	m4_esyscmd_s([
		git for-each-ref --sort='-*committerdate' --format \
			'News for %(refname:short):%0a===============%0a%0a%(body)' \
			refs/tags > NEWS
	])
	, []
)

# read version number
define([gitversion], m4_esyscmd_s([cat VERSION]))

# initialize configure
AC_INIT([hep-mc], [gitversion])

# initialize automake
AM_INIT_AUTOMAKE([-Wall -Werror dist-bzip2 no-dist-gzip nostdinc 1.11])

# use silent rules by default
AM_SILENT_RULES([yes])

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

AC_ARG_ENABLE([examples], AS_HELP_STRING([--enable-examples],
	[Enable building of example programs]))

AM_CONDITIONAL([HAVE_EXAMPLES], [test "x$enable_examples" = "xyes"])

AS_IF([test "x$enable_examples" = "xyes"], [
	LX_FIND_MPI
])

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
