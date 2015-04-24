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
- *compatible* to the new C++11 random number generators interface, see e.g.
  here: http://en.cppreference.com/w/cpp/numeric/random. This enables one to
  choose between a variety of standard generators or to write a custom one and
  use it with the integration routines.

Furthermore, for VEGAS

- one can specify the sample size for each iteration separately; afterwards
- one can selectively include/exclude each iteration's result in a cumulative
  estimate,
- it is possible to set the bin count to study its effect on the precision for
  complicated integrands, for example. In other VEGAS implementations the bin
  count is set to a fixed number, e.g. 50 in `GSL`_ and 128 in `CUBA`_.

Installation
============

If you downloaded the sources with ``git`` or as a ZIP-archive you first have to
generate the autotools files. Run::

    autoreconf -fiv

inside the project's directory. This will generate the missing files, e.g. the
configure script and makefiles. Alternatively, you can download a pre-built
tarball from the `releases page <http://github.com/cschwan/hep-mc/releases>`_.

To install it, use the usual sequence of commands for autotools-based projects,
for example::

    ./configure --prefix=/usr/local
    make
    make install

The ``prefix`` specifies the directory in which the library will be installed
to. Make sure your compiler will find this directory.

In addition you may pass the ``--enable-doxygen`` option to ``./configure`` to
build the HTML documentation. For this to work you will need a recent version of
the `Doxygen`_ program.

A complete list of configuration options is available by typing::

    ./configure --help

Documentation
=============

Documentation is available online at http://cschwan.github.io/hep-mc or can be
generated from sources with Doxygen (see Installation_).

Usage
=====

Since this library uses features from the new C++11 standard, you have to enable
these with your compiler. For GCC and clang this can be done by passing an
additional parameter to the compiler, e.g.

    g++ -std=c++0x my_program.cpp

Since this project is based on templates, the inclusion of the main header,

    #include <hep/mc.hpp>

is sufficient to use it and you do not need to link against a library. To see
the library in action take a look the example programs in the
`examples directory`_.

.. _GSL: http://www.gnu.org/software/gsl/
.. _CUBA: http://www.feynarts.de/cuba/
.. _Doxygen: http://www.doxygen.org/
.. _examples directory: http://github.com/cschwan/hep-mc/tree/master/examples
