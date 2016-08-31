GENERAL
=======

- add ``non_zero_calls()`` member function to ``_result`` classes
- add serialization functions for the ``_result`` classes that are needed to
  exchange and accumulate the results from different MPI processors
- write cumulative functions that accumulate also the distributions
- remove the globally set callback functions and make them local to every call
  to an integrator?
- add functionality to capture non-finite numbers of the integrand and report
  how many numbers of the generator need to be discarded to get there: save a
  copy of the generator the integrator obtains, print out the state (with
  ``operator<<``) and the number of discards needed
- save largest (and smallest?) weight?

NEW INTEGRATORS
===============

- write new integrator using the FOAM algorithm?

MULTI CHANNEL
=============

- add possibility for the controlled removal of channels?
- add beta parameter to dampen the weight adjustment for better convergence

VEGAS
=====

- get rid of the global configuration option, make it an argument to the vegas
  function?
