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
	$(INSTALL) -d $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) ppidd-config $(DESTDIR)$(bindir)
	$(INSTALL) -d $(DESTDIR)$(libdir)
	$(INSTALL_DATA) libppidd.a $(DESTDIR)$(libdir)
	$(INSTALL) -d $(DESTDIR)$(includedir)
	$(INSTALL_DATA) src/ppidd.h $(DESTDIR)$(includedir)
	$(INSTALL_DATA) src/mpimutex.h $(DESTDIR)$(includedir)
	$(INSTALL_DATA) src/mpiga_base.h $(DESTDIR)$(includedir)
ifdef FC
	$(INSTALL_DATA) src/$(subst %M,PPIDD,$(subst %m,ppidd,$(FC_MODFMT))) $(DESTDIR)$(includedir)
endif

.PHONY: clean
clean:
	git clean -X -d -f
