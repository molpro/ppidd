include config.mk

.PHONY: default
default: libppidd.a

libppidd.a:
	@$(MAKE) -C src
	@rm -f $@
	@$(AR) $(ARFLAGS) $@ src/*.o
	@$(RANLIB) $@

.PHONY: test
test: library
	@echo
	@echo 'Building PPIDD test suite'
	@echo
	$(MAKE) -C test

.PHONY: doc
doc:
ifdef DOXYGEN
	@cd doc; $(DOXYGEN) Doxyfile
else
	@echo 'doxygen does not appear to be installed'
endif

clean:
	@$(foreach directory,src test,$(MAKE) -C $(directory) clean;)
	@rm -rf libppidd.a
