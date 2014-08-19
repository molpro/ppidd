include config.mk

#ifeq ($(notdir $(FC)),ifort)
#override FFLAGS+=-Vaxlib
#endif

#LIBS+=-L$(realpath $(wildcard $(firstword $(INCLUDE))/../lib)) -lga -larmci

.PHONY: default
default: library

.PHONY: library
library:
	@rm -rf lib
	@mkdir lib
	$(MAKE) -C src
	@$(AR) $(ARFLAGS) lib/libppidd.a src/*.o
	@$(RANLIB) lib/libppidd.a

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
	@rm -rf lib
