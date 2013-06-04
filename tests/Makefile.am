AUTOMAKE_OPTIONS = parallel-tests

AM_CPPFLAGS = \
	-I$(top_srcdir)/include \
	-DBOOST_TEST_NO_LIB \
	-DBOOST_TEST_DYN_LINK \
	-DBOOST_TEST_MAIN

AM_DEFAULT_SOURCE_EXT = .cpp

check_PROGRAMS = \
	test_plain \
	test_vegas

LDADD = -lboost_unit_test_framework

TESTS = $(check_PROGRAMS)
