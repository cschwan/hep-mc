html.tar: Doxyfile references.bib $(top_builddir)/include/stamp-headers
	$(AM_V_GEN)$(DOXYGEN) Doxyfile
	$(AM_V_at)tar cf html.tar html/

all-local: html.tar

mostlyclean-local:
	rm -rf html/

install-data-local:
	$(MKDIR_P) "$(DESTDIR)$(docdir)"
	tar xf html.tar -C "$(DESTDIR)$(docdir)"

uninstall-local:
	for i in `tar tf html.tar`; do \
	    if [ -f "$(DESTDIR)$(docdir)"/$$i ]; then \
	        rm "$(DESTDIR)$(docdir)"/$$i; \
	    fi \
	done

EXTRA_DIST = references.bib
CLEANFILES = html.tar
