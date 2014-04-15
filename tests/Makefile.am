AM_CPPFLAGS = -I$(top_srcdir)/include -I$(srcdir)/3rd-party

EXTRA_DIST = 3rd-party/catch.hpp

AM_DEFAULT_SOURCE_EXT = .cpp

check_PROGRAMS = \
	test_plain \
	test_vegas

# redundant for automake >= 1.11
test_plain_SOURCES = test_plain.cpp
test_vegas_SOURCES = test_vegas.cpp

TESTS = $(check_PROGRAMS)
