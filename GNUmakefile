include config.mk

.PHONY: default
default: libppidd.a

libppidd.a: $(wildcard src/*.cpp src/*.h src/*.F) ppiddcxx ppiddfc
	@$(MAKE) -C src
	@rm -f $@
	@$(AR) $(ARFLAGS) $@ src/*.o
	@$(RANLIB) $@

ppiddcxx: config.mk
	@echo '#/bin/sh' > $@
	@echo '$(CXX) $(CPPFLAGS) $$*' >> $@
	@chmod +x $@

ifeq ($(FC),)
.PHONY: ppiddfc
ppiddfc:
else
ppiddfc: config.mk
	@echo '#/bin/sh' > $@
	@echo '$(FC) $(CPPFLAGS) $$*' >> $@
	@chmod +x $@
endif

.PHONY: test
test: libppidd.a
	@echo
	@echo 'Building PPIDD test suite'
	@echo
	@$(MAKE) -C test

.PHONY: doc
doc:
ifdef DOXYGEN
	@$(DOXYGEN) src/Doxyfile 1>$@.log
else
	@echo 'doxygen does not appear to be installed'
endif

.PHONY: install
install: default
	@mkdir -p $(prefix)/bin
	@cp -v ppiddcxx $(prefix)/bin/
ifneq ($(FC),)
	@cp -v ppiddfc $(prefix)/bin/
endif
	@mkdir -p $(prefix)/lib
	@cp -v libppidd.a $(prefix)/lib/

.PHONY: clean
clean:
	git clean -X -d -f
