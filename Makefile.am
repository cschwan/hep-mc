ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = NEWS.rst README.rst

# if doxygen is enabled include doc subdir
if HAVE_DOXYGEN
    DOCDIR = doc
endif

# if examples are enabled include examples subdir
if HAVE_EXAMPLES
    EXAMPLESDIR = examples
endif

SUBDIRS = include $(DOCDIR) $(EXAMPLESDIR) tests
