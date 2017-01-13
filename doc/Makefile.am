docpkg = $(PACKAGE_TARNAME)-doc-$(PACKAGE_VERSION).tar

CLEANFILES = $(docpkg) doxygen_sqlite3.db

EXTRA_DIST = \
	callbacks.dox \
	distributions.dox \
	Doxyfile.in \
	DoxygenLayout.xml \
	examples.dox \
	extra.css \
	footer.html \
	global_configuration.dox \
	integrands.dox \
	mainpage.dox \
	multi_channel.dox \
	plain.dox \
	references.bib \
	results.dox \
	vegas.dox \
	vegas_pdf.dox

$(docpkg): Doxyfile $(EXTRA_DIST) $(top_builddir)/include/stamp-headers
	$(AM_V_GEN)$(DOXYGEN) Doxyfile
	$(AM_V_at)tar cf $(docpkg) html/

all-local: $(docpkg)

mostlyclean-local:
	rm -rf html/

install-data-local:
	$(MKDIR_P) "$(DESTDIR)$(docdir)"
	tar xf $(docpkg) -C "$(DESTDIR)$(docdir)"

uninstall-local:
	for i in `tar tf $(docpkg)`; do \
	    if [ -f "$(DESTDIR)$(docdir)"/$$i ]; then \
	        rm "$(DESTDIR)$(docdir)"/$$i; \
	    fi \
	done
