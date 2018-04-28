include ../config.mk

OBJECTS=$(addsuffix .o,$(basename $(wildcard *.cpp)))
ifneq ($(FC),)
OBJECTS+=ppidd_module.o
endif

.PHONY: default
default: $(OBJECTS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) -c $(CXXFLAGS) $<

%.o: %.F
	$(FC) -c $(FCFLAGS) $<
