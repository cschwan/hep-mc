New in 0.2.2.9999:
==================

- add more unit tests, use Google Test as testing framework
- improved the MPI VEGAS example
- added new ``hep::mpi_single_generator`` function that adds the possibility to
  use the same random numbers for MPI and non-MPI routines. This makes
  ``hep::mpi_vegas`` return the same numerical result as ``hep::vegas``
- fixed integer overflow bug on 32-bit architectures that caused VEGAS' grid
  adaption to fails when the number of calls surpassed 2^16
- modified VEGAS callback functions to return a boolean variable signaling to
  stop all remaining iterations
- modified the verbose callback function to print the error in percentage

New in 0.2.2:
=============

- added two examples showing the advanced use of the VEGAS integration routine
- added a callback mechanism for the VEGAS routines; this enables one to print
  intermediate results from completed iterations
- added VEGAS routines that reuse existing grids
- grids can now be easily saved to files with iostream shift operators, see the
  example ``examples/read_linear_grid.cpp`` which can be used with the grid
  file in ``examples/grid_file``


New in 0.2.1:
=============

- added MPI-parallelized PLAIN integration algorithm
- lowered dependency on autotools features

New in 0.2:
===========

- added new functions for combining several VEGAS results
- added MPI-based paralellized VEGAS integration algorithm
- improved documentation
- improved VEGAS performance

New in 0.1:
===========

- added PLAIN integrator
- added VEGAS integrator
