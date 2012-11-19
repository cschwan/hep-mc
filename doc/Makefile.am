doxygen-files: Doxyfile references.bib $(top_builddir)/include/stamp-headers
	$(DOXYGEN) Doxyfile
	tar cf html.tar html/
	tar tf html.tar > doxygen-files

all-local: doxygen-files

mostlyclean-local:
	rm -rf html/

install-data-local:
	$(MKDIR_P) "$(DESTDIR)$(docdir)"
	cd "$(DESTDIR)$(docdir)"; tar xof "$(abs_builddir)"/html.tar

uninstall-local:
	for i in `cat doxygen-files`; do \
	    if [ -f "$(DESTDIR)$(docdir)"/$$i ]; then \
	        rm "$(DESTDIR)$(docdir)"/$$i; \
	    fi \
	done

EXTRA_DIST = references.bib
MOSTLYCLEANFILES = doxygen-files
CLEANFILES = html.tar
