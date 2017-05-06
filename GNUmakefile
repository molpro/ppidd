include config.mk

.PHONY: default
default: libppidd.a

libppidd.a: $(wildcard src/*.cpp src/*.h src/*.F)
	@$(MAKE) -C src
	@rm -f $@
	@$(AR) $(ARFLAGS) $@ src/*.o
	@$(RANLIB) $@

.PHONY: test
test: libppidd.a
	@$(MAKE) -C $@

.PHONY: doc
doc:
ifdef DOXYGEN
	@$(DOXYGEN) src/Doxyfile 1>$@.log
else
	@echo 'doxygen does not appear to be installed'
endif

.PHONY: install
install: default
	$(INSTALL) -d $(DESTDIR)$(libdir)
	$(INSTALL_DATA) -t $(DESTDIR)$(libdir) libppidd.a
	$(INSTALL) -d $(DESTDIR)$(includedir)
	$(INSTALL_DATA) -t $(DESTDIR)$(includedir) src/ppidd_c.h
	$(INSTALL_DATA) -t $(DESTDIR)$(includedir) src/ppidd_eaf_c.h
	$(INSTALL_DATA) -t $(DESTDIR)$(includedir) src/ppidd_sf_c.h

.PHONY: clean
clean:
	git clean -X -d -f
