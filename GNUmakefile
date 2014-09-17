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

clean:
	@$(foreach directory,src test,$(MAKE) -C $(directory) clean;)
	@rm -rf libppidd.a
