Project Description
===================

``hep-mc`` is a C++11 template library providing Monte Carlo integration
algorithms, currently only PLAIN and VEGAS. In addition, this library comes with
MPI-parallelized routines.

The aim of this project is to provide functions and classes that are:

- *easily usable*, e.g. the default VEGAS routine can be called with only three
  parameters. More can be given if the default parameters are not satisfactory,
- *templatized* in order to support a wide range of floating point types used to
  perform the numerical computations (e.g. ``float``, ``double``,
  ``long double`` and custom types),
- *modularized* which makes it easy to develop new integration algorithms on top
  of known ones - e.g. VEGAS with a modified rebinning algorithm,
- *compatible* to the new C++11 random number generators [1]_. As a side effect,
  it is possible to write and use custom generators, e.g. to generate
  quasi-random numbers.

Furthermore, for VEGAS

- one can specify the sample size for each iteration separately; afterwards
- one can selectively include/exclude each iteration's result in a cumulative
  estimate,
- it is possible to set the bin count to study its effect on the precision for
  complicated integrands, for example. In other VEGAS implementations the bin
  count is set to a fixed number, e.g. 50 in GSL and 128 in CUBA.

Installation and Testing
========================

``hep-mc`` is a header-only library which requires no compilation. To install
it, use the usual sequence of commands for autotools-based projects, for
example::

    ./configure --prefix=/usr/local --enable-doxygen
    make
    make install

The ``prefix`` specifies the directory in which the library will be installed
to. If you enable Doxygen, ``make`` will generate the HTML documentation for
this library which will be installed as well. Note that you need a recent
version of this program [2]_.

A complete list of configuration options is available by typing::

    ./configure --help

Missing ``./configure``
=======================

If you have obtained this repository which does not contain the files generated
by autotools ...

- run ``autoreconf -fiv`` to generate the files which are not under revision
  control
- if the last step fails, your autotools are most likely outdated. Update
  autoconf and automake and try again.

.. [1] See e.g. http://en.cppreference.com/w/cpp/numeric/random
.. [2] http://www.doxygen.org/
