AUTOMAKE_OPTIONS = parallel-tests

AM_CPPFLAGS = \
	-I$(top_srcdir)/include \
	-DBOOST_TEST_NO_LIB \
	-DBOOST_TEST_DYN_LINK \
	-DBOOST_TEST_MAIN

AM_DEFAULT_SOURCE_EXT = .cpp

check_PROGRAMS = test_plain

LDADD = -lboost_unit_test_framework

LOG_COMPILER = BOOST_TEST_BUILD_INFO=yes BOOST_TEST_DETECT_FP_EXCEPTIONS=yes

TESTS = $(check_PROGRAMS)
