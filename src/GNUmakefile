include ../config.mk

ifeq ($(FC),)
SOURCEFILES+=$(filter-out %fortran.cpp,$(wildcard *.cpp))
else
SOURCEFILES+=$(wildcard *.cpp *.F)
endif

.PHONY: default
default: $(addsuffix .o,$(basename $(SOURCEFILES)))

ppidd_fortran.o ppidd_c.o: ppidd_share.h

%.o: %.cpp
	$(CXX) $(CPPFLAGS) -c $(CXXFLAGS) $<

%.o: %.F
	$(FC) -c $(FFLAGS) $<
