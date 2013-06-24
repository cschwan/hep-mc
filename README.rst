Project Description
===================

``hep-mc`` is a C++11 template library providing Monte Carlo integration
algorithms, currently only PLAIN and VEGAS.

Installation and Testing
========================

hep-mc is a header-only library which requires no compilation. To install it,
use the usual sequence of commands for autotools-based projects, for example::

    ./configure --prefix=/usr/local --enable-doxygen
    make
    make install

This will generate the Doxygen documentation for this library first and installs
it in the second step to /usr/local/doc/hep-ga. The headers are installed to
/usr/local/include. You may change these directories by specifying additional
options to ./configure. The following command gives you a complete list of
available configuration options::

    ./configure --help

Missing ``./configure``
=======================

If you have obtained this repository which does not contain the files generated
by autotools ...

- run ``autoreconf -fiv`` to generate the files which are not under revision
  control
- if the last step failed, your autotools are most likely outdated. Update
  autoconf and automake and try again.

Usage
=====

See Doxygen documentation.
