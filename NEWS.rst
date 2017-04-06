New in 0.5:
===========

- fixed a bug that caused NaNs to slip through to the integrated result if they
  were coming from the point's weight
- it is now possible to disable channels by setting their weight to zero
- the function ``multi_channel_verbose_callback`` now prints an overview of the
  smallest and largest a-priori weights
- the functions ``multi_channel`` and ``mpi_multi_channel`` now accept an
  optional parameter ``min_calls_per_channel`` that limits the a-priori weights
  to a smallest weight determined by a minimum number of calls. If set to zero
  (the default), this corresponds to the old behavior
- WARNING: this version introduces many interface changes; if you used previous
  version of this library you must update your code accordingly
- implemented lazy evaluation of the weight for the multi channel integrator.
  This required to change the interface of the map function which now must
  accept an additional parameter ``action`` that determines whether the
  ``coordinates`` or the ``densities`` must be calculated. This in turn enables
  the algorithm to skip the evaluation of ``densities`` when the integrand
  returns a zero value
- changed name of ``density_function()`` in ``hep::multi_channel_point`` to
  ``map()``
- when calling the integration algorithms the numeric type must no longer be
  explicitly specified. The type information is now infered from the integrand
- renamed the class ``vegas_iteration_result`` to ``vegas_result``
- completely changed the interfaces of all integrators; the functions that are
  integrated must now be specified by using an ``hep::integrand`` or
  ``hep::multi_channel_integrand`` which are obtained by using
  ``hep::make_integrand`` or ``hep::make_multi_channel_integrand``. This allows
  the library to support the generation of differential distributions for all
  integrators. If distributions should be generated, the integrands must accept
  a second parameter, a reference to ``hep::projector``, that takes care of the
  binning. The results of the distributions are captured in a
  ``hep::distribution_result`` for each distribution
- ``hep::plain`` now returns a ``hep::plain_result`` which contains the results
  for possibly generated distributions
- all ``_result`` classes except ``hep::mc_result`` derive now from
  ``hep::plain_result``
- encapsulated all members of the class ``mc_point`` and classes deriving from
  it
- WARNING: interface changes to the class ``mc_point`` and all deriving classes;
  the member function does no longer include the averaging factor 1/N where N is
  the number of calls. If you used this member variable before, you MUST CHANGE
  YOUR CODE accordingly

New in 0.4:
===========

- added multi channel MC integrator with example
- added a new MC helper function ``hep::multi_channel_max_difference`` that
  calculates the maximum difference for adaptive optimization in multi channel
  MC
- added callbacks for multi channel MC analogous to the case for VEGAS
- added MPI parallelized routines for multi channel MC

New in 0.3.1:
=============

- added new internal class ``hep::kahan_accumulator``. This class takes care of
  the kahan summation
- changed interfaces ``mc_result`` and ``vegas_iteration_result``: Their member
  variables are now properly encapsulated and accessible via public member
  functions. Furthermore, the ``mc_result`` now allows to obtain ``sum()``,
  ``sum_of_squares()``, and ``variance()``

New in 0.3:
===========

- fixed compilation error
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
