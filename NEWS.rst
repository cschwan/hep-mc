New in 0.2.1.9999:
==================

- added two examples showing the advanced use of the VEGAS integration routine
- added a callback mechanism for the VEGAS routines; this enables one to print
  intermediate results from completed iterations
- added VEGAS routines to reuse an existing grid
- grids can now be easily saved to files with iostream shift operators, see the
  example 'examples/read_linear_grid.cpp' which can be used with the grid file
  in 'examples/grid_file'

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
