AM_CPPFLAGS = \
	-I$(top_srcdir)/include \
	-DBOOST_TEST_NO_LIB \
	-DBOOST_TEST_DYN_LINK \
	-DBOOST_TEST_MAIN

AM_DEFAULT_SOURCE_EXT = .cpp

check_PROGRAMS = \
	test_plain \
	test_vegas

# redundant for automake >= 1.11
test_plain_SOURCES = test_plain.cpp
test_vegas_SOURCES = test_vegas.cpp

LDADD = -lboost_unit_test_framework

TESTS = $(check_PROGRAMS)
