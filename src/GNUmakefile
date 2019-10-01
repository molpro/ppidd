include ../config.mk
VPATH=$(ppidd_srcdir)/src

OBJECTS=$(notdir $(addsuffix .o,$(basename $(wildcard $(VPATH)/*.cpp))))
ifneq ($(FC),)
OBJECTS+=ppidd_module.o
endif

.PHONY: default
default: $(OBJECTS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) -I. -c $(CXXFLAGS) $<

%.o: %.F90
	$(FC) -c $(FCFLAGS) $<
