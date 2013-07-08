AM_CPPFLAGS = -I$(top_srcdir)/include

AM_DEFAULT_SOURCE_EXT = .cpp

noinst_PROGRAMS = \
	read_linear_grid \
	mpi_vegas_example \
	vegas_example \
	vegas_grid

# redundant for automake >= 1.11
read_linear_grid_SOURCES = read_linear_grid.cpp
mpi_vegas_example_SOURCES = mpi_vegas_example.cpp
vegas_example_SOURCES = vegas_example.cpp
vegas_grid_SOURCES = vegas_grid.cpp

mpi_vegas_example_CXXFLAGS = $(MPI_CXXFLAGS)
mpi_vegas_example_LDADD = $(MPI_CXXLDFLAGS)
