AM_CPPFLAGS = -I$(top_srcdir)/include -I$(srcdir) -DGTEST_HAS_PTHREAD=0

AM_DEFAULT_SOURCE_EXT = .cpp

# convenience libraries for the google unit test framework
noinst_LTLIBRARIES = libgtest.la libgtest_main.la
libgtest_la_SOURCES = gtest/gtest.h gtest-all.cpp
libgtest_main_la_SOURCES = gtest/gtest.h gtest_main.cpp

check_PROGRAMS = \
	test_mc_helper \
	test_mc_point \
	test_mc_result \
	test_piecewise_constant_pdf \
	test_plain \
	test_vegas

# redundant for automake >= 1.11
test_mc_helper_SOURCES = test_mc_helper.cpp
test_mc_point_SOURCES = test_mc_point.cpp
test_mc_result_SOURCES = test_mc_result.cpp
test_piecewise_constant_pdf_SOURCES = test_piecewise_constant_pdf.cpp
test_plain_SOURCES = test_plain.cpp
test_vegas_SOURCES = test_vegas.cpp

LDADD = $(noinst_LTLIBRARIES)
TESTS = $(check_PROGRAMS)
