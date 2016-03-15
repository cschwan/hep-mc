GENERAL
=======

- make callback functions argument to the algorithms?
- write cumulative functions that accumulate also the distributions

MULTI CHANNEL
=============

- set a minimal bound on the a-priori weights?

  + The number of calls for each channel must not be equal to one, because a
    single call means an infinitely large error.
  + If the number of calls is exactly zero (which can easily happen if either
	the weight crosses a minimal threshold or if the number of channels is too
	large in comparison with the total number of call), the channel is
    effectively removed.
  + If the number of calls is small but larger than one the channel will mainly
	contribute error to the integration
  + This leads to the question: What is the smallest number of calls for each
    channel, or equivalently, what is the minimal weight for each channel?

- add possibility for the controlled removal of channels?
- add beta parameter to dampen the weight adjustment for better convergence
- rename ``multi_channel_point2::density_function()`` to something shorther,
  ``map()`` maybe?

VEGAS
=====

- get rid of the global configuration option, make it aan argument to the vegas
  function?
