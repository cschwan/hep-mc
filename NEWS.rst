New in 0.3:
===========

- removal of the ``hep::mpi_single_generator`` function. This is now on by
  default - the integration results are independent of the number of processes
- temporary fix for single-precision floating point numbers that are ``1.0f``
  and therefore outside of the half-open interval [0,1)
- changed ``hep::vegas_pdf`` interface to use more descriptive member function
  names and moved the member function ``icdf`` outside the class
- revamped documentation
- changed function names of ``hep::cumulative_result`` and
  ``hep::chi_square_dof`` to the same names with a ``0`` and ``1`` at the end
  of their names; ``0`` uses the same algorithm as before, ``1`` weighs all
  results equally. The function ``hep::cumulative_result`` remains but requires
  an additional argument
- added new global configuration function ``hep::vegas_cuba_refinement()``
  that can be used to modify pdf refinement to use CUBA's method
- removed ``--enable-tests`` configure switch. All external dependecies are
  gone and running ``make check`` is sufficient now to run tests
- changed interfaces ``hep::cumulative_result`` and ``hep::chi_square_dof``
  which no longer need to be called with the numeric type; the type is
  automatically determined by the iterators
- support for MPI can now be enabled with ``--enable-mpi`` independently from
  the examples. If MPI is enabled, the MPI examples are built if examples are
  activated and the MPI tests are checked if tests are activated. The MPI
  headers are installed in any case because MPI must be explicity requested by
  including ``hep/mc-mpi.hpp``
- renamed ``hep::linear_grid`` to ``hep::vegas_pdf`` and moved code from VEGAS
  into this class
- add more unit tests, use Google Test as testing framework
- improved the MPI VEGAS example
- added new ``hep::mpi_single_generator`` function that adds the possibility to
  use the same random numbers for MPI and non-MPI routines. This makes
  ``hep::mpi_vegas`` return the same numerical result as ``hep::vegas``
- fixed integer overflow bug on 32-bit architectures that caused VEGAS' grid
  adjustment to fail when the number of calls surpassed 2^16
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
