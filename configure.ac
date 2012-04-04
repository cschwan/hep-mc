AC_PREREQ([2.65])

# check if git is available and current directory contains repository
m4_esyscmd_s([git status])

# if repository is available use git to create files which must be available
# before AC_INIT is called:
ifelse(sysval,
	# VERSION file, ...
	m4_esyscmd_s([echo `git describe 2>/dev/null` | tee VERSION])

	# AUTHORS file and ...
	m4_esyscmd_s([git log --format='%aN <%aE>' | uniq > AUTHORS])

	# ChangeLog file
	m4_esyscmd_s(
		[git log --date=short --format='%ad  %aN  <%aE> %n%x09* %s%n' \
		> ChangeLog]
	),
	[]
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

# point to a source file in order to let ./configure find itself
#AC_CONFIG_SRCDIR([])

# do not pollute command line with compiler definitions but instead write them
# into the specified header
AC_CONFIG_HEADERS([config.hpp])

AC_ARG_ENABLE(
	[examples],
	AS_HELP_STRING([--enable-examples],
		[Enable building and installation (only sources) of example programs])
)

AM_CONDITIONAL([HAVE_EXAMPLES], [test "x$enable_examples" = "xyes"])

# add possibility to generate API documentation with Doxygen
AC_ARG_ENABLE(
	[doxygen],
	AS_HELP_STRING([--enable-doxygen],
		[Enable generation of Doxygen documentation])
)

# variable pointing to the doxygen binary
AC_ARG_VAR([DOXYGEN], [Documentation system])

# check for the doxygen binary
AS_IF([test "x$enable_doxygen" = "xyes"], [
	AC_CHECK_PROGS([DOXYGEN], [doxygen])
	AS_IF(
		[test -z "$DOXYGEN"],
		AC_MSG_ERROR([Doxygen not found but requested])
	)
])

# if doxygen is requested, define HAVE_DOXYGEN and configure Doxyfile.in
AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])

# check for boost test framework library
AC_SEARCH_LIBS(
	[_ZN5boost9unit_test14unit_test_mainEPFbvEiPPc],
	[boost_unit_test_framework],
	[],
	AC_MSG_WARN([
    --------------------------------------------------
    Did not find Boosts Unit Test Framework library.
    You will not be able to compile and run the tests.
    --------------------------------------------------])
)

# check for the corresponding header
AC_CHECK_HEADERS(
	[boost/test/unit_test.hpp],
	[],
	AC_MSG_WARN([
    --------------------------------------------------
    Did not find Boosts Unit Test Framework header.
    You will not be able to compile and run the tests.
    --------------------------------------------------])
)

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
