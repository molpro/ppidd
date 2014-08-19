include ../../config.mk
AUTOMAKE_FILES=config.guess config.sub install-sh

.PHONY: default
default: ../../configure $(AUTOMAKE_FILES)

../../configure: configure.ac $(wildcard *.m4)
ifeq ($(AUTOCONF),)
	@echo "AUTOCONF undefined, $@ not rebuilt"
else
	@$(AUTOCONF) -Wall $< > $@
	@chmod +x $@
	@rm -rf autom4te.cache
endif

.PHONY: $(AUTOMAKE_FILES)
$(AUTOMAKE_FILES):
ifeq ($(AUTOMAKE_LIBDIR),)
	@echo "AUTOMAKE_LIBDIR undefined, $@ not rebuilt"
else
	@cp -v $(AUTOMAKE_LIBDIR)/$@ .
endif
