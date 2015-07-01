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

Usage
=====

Since this library uses features from the new C++11 standard, you have to enable
these with your compiler. For the GCC and clang compilers this can be done by
passing an additional parameter to the compiler, e.g.::

    g++ -std=c++0x my_program.cpp

Since this project is based on templates, the inclusion of the main header,::

    #include <hep/mc.hpp>

is sufficient to use it and you do not need to link against a library. To see
the library in action take a look the example programs in the
`examples directory`_.

Documentation and Examples
==========================

Documentation is available online at http://cschwan.github.io/hep-mc or can be
generated from sources (see Installation_). The examples can be viewed from
within the documentation.

Installation
============

The easiest way to use this library is to just download it from the `releases
page`_ and point your compiler to the `include directory`_ - there is no library
that needs to be compiled.

If you want to automatically compile the example programs, generate the
documentation, and/or install the headers you can also use the usual sequence
for autotools-based projects, i.e.::

    ./configure --prefix=/usr/local
    make
    make install

which installs the headers into ``/usr/local/include``. To generate the
documentation pass ``--enable-doxygen`` to ``./configure``. Example programs are
compiled when you add ``--enable-examples``. For more options type ``./configure
--help``.

When you got the sources with ``git`` or downloaded them as a ZIP file, you will
notice that the configure script is missing. Then run::

    autoreconf -fiv

inside the project's top-level directory. This will generate the missing files.

.. _releases page: http://github.com/cschwan/hep-mc/releases
.. _include directory: http://github.com/cschwan/hep-mc/tree/master/include
.. _examples directory: http://github.com/cschwan/hep-mc/tree/master/examples
