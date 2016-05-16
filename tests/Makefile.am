AM_CPPFLAGS = -I$(top_srcdir)/include -I$(srcdir) -DGTEST_HAS_PTHREAD=0

AM_DEFAULT_SOURCE_EXT = .cpp

# convenience libraries
check_LIBRARIES = libgenz.a libgtest.a libgtest_main.a

libgenz_a_SOURCES = genz_integrand.hpp genz_integrand.cpp
libgtest_a_SOURCES = gtest/gtest.h gtest-all.cpp
libgtest_main_a_SOURCES = gtest/gtest.h gtest_main.cpp

serial_tests = \
	test_discrete_distribution \
	test_distributions \
	test_genz_integrand \
	test_mc_helper \
	test_mc_point \
	test_mc_result \
	test_numerical_results \
	test_vegas_icdf_bug \
	test_vegas_pdf

check_PROGRAMS = $(serial_tests)

LDADD = libgtest.a libgtest_main.a

test_genz_integrand_LDADD = libgenz.a $(LDADD)

TESTS = $(serial_tests)

if HAVE_MPI

check_LIBRARIES += libgtest_mpi_main.a

libgtest_mpi_main_a_SOURCES = gtest/gtest.h gtest_mpi_main.cpp
libgtest_mpi_main_a_CPPFLAGS = $(AM_CPPFLAGS) $(MPI_CPPFLAGS)

parallel_scripts = \
	test_mpi_numerical_results.sh

parallel_tests = \
	test_mpi_numerical_results

check_PROGRAMS += $(parallel_tests)
nodist_check_SCRIPTS = $(parallel_scripts)
TESTS += $(parallel_scripts)

test_mpi_numerical_results_CPPFLAGS = \
	$(AM_CPPFLAGS) $(MPI_CPPFLAGS) -DHEP_USE_MPI
test_mpi_numerical_results_LDFLAGS = $(MPI_LDFLAGS)
test_mpi_numerical_results_LDADD = libgtest.a libgtest_mpi_main.a $(MPI_LIBS)
test_mpi_numerical_results_SOURCES = test_numerical_results.cpp

endif

