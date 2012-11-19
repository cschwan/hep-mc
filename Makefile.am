ACLOCAL_AMFLAGS = -I m4 --install

EXTRA_DIST = VERSION

# if doxygen is enabled include doc subdir
if HAVE_DOXYGEN
    DOCDIR = doc
endif

# if examples are enabled include examples subdir
if HAVE_EXAMPLES
    EXAMPLESDIR = examples
endif

# if tests are enabled include tests subdir
if HAVE_TESTS
    TESTSDIR = tests
endif

SUBDIRS = include $(DOCDIR) $(EXAMPLESDIR) $(TESTSDIR)
