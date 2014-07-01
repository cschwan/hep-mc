AM_CPPFLAGS = -I$(top_srcdir)/include -I$(srcdir) -DGTEST_HAS_PTHREAD=0

AM_DEFAULT_SOURCE_EXT = .cpp

# convenience libraries for the google unit test framework
check_LIBRARIES = libgtest.a libgtest_main.a
libgtest_a_SOURCES = gtest/gtest.h gtest-all.cpp
libgtest_main_a_SOURCES = gtest/gtest.h gtest_main.cpp

check_PROGRAMS = \
	test_mc_helper \
	test_mc_point \
	test_mc_result \
	test_plain \
	test_vegas \
	test_vegas_pdf

# redundant for automake >= 1.11
test_mc_helper_SOURCES = test_mc_helper.cpp
test_mc_point_SOURCES = test_mc_point.cpp
test_mc_result_SOURCES = test_mc_result.cpp
test_plain_SOURCES = test_plain.cpp
test_vegas_SOURCES = test_vegas.cpp
test_vegas_pdf_SOURCES = test_vegas_pdf.cpp

LDADD = libgtest.a libgtest_main.a

if HAVE_MPI

check_LIBRARIES += libgtest_mpi_main.a
libgtest_mpi_main_a_SOURCES = gtest/gtest.h gtest_mpi_main.cpp
libgtest_mpi_main_a_CPPFLAGS = $(AM_CPPFLAGS) $(MPI_CPPFLAGS)

check_PROGRAMS += \
	test_mpi_vegas

test_mpi_vegas_CPPFLAGS = $(AM_CPPFLAGS) $(MPI_CPPFLAGS)
test_mpi_vegas_LDFLAGS = $(MPI_LDFLAGS)
test_mpi_vegas_LDADD = libgtest.a libgtest_mpi_main.a $(MPI_LIBS)

# redundant for automake >= 1.11
test_mpi_vegas_SOURCES = test_mpi_vegas.cpp

endif

TESTS = $(check_PROGRAMS)
