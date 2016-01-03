AM_CPPFLAGS = -I$(top_srcdir)/include

AM_DEFAULT_SOURCE_EXT = .cpp

noinst_PROGRAMS = \
	distributions_example \
	helper_example \
	multi_channel_example \
	vegas_example \
	vegas_read_pdf \
	vegas_write_pdf \
	vegas_stop_after_precision

if HAVE_MPI

noinst_PROGRAMS += \
	mpi_vegas_example \
	mpi_multi_channel_example

mpi_vegas_example_CPPFLAGS = $(AM_CPPFLAGS) $(MPI_CPPFLAGS)
mpi_vegas_example_LDADD = $(MPI_LIBS)
mpi_vegas_example_LDFLAGS = $(MPI_LDFLAGS)

mpi_multi_channel_example_CPPFLAGS = $(AM_CPPFLAGS) $(MPI_CPPFLAGS) -DUSE_MPI
mpi_multi_channel_example_LDADD = $(MPI_LIBS)
mpi_multi_channel_example_LDFLAGS = $(MPI_LDFLAGS)
mpi_multi_channel_example_SOURCES = multi_channel_example.cpp

# redundant for automake >= 1.11
mpi_vegas_example_SOURCES = mpi_vegas_example.cpp

endif
