AM_CPPFLAGS = -I$(top_srcdir)/include

AM_DEFAULT_SOURCE_EXT = .cpp

noinst_PROGRAMS = \
	helper_example \
	read_linear_grid \
	vegas_example \
	vegas_grid \
	vegas_stop_after_precision

# redundant for automake >= 1.11
helper_example_SOURCES = helper_example.cpp
read_linear_grid_SOURCES = read_linear_grid.cpp
vegas_example_SOURCES = vegas_example.cpp
vegas_grid_SOURCES = vegas_grid.cpp
vegas_stop_after_precision_SOURCES = vegas_stop_after_precision.cpp

if HAVE_MPI

noinst_PROGRAMS += \
	mpi_vegas_example

mpi_vegas_example_CPPFLAGS = $(AM_CPPFLAGS) $(MPI_CPPFLAGS)
mpi_vegas_example_LDADD = $(MPI_LIBS)
mpi_vegas_example_LDFLAGS = $(MPI_LDFLAGS)

# redundant for automake >= 1.11
mpi_vegas_example_SOURCES = mpi_vegas_example.cpp

endif
