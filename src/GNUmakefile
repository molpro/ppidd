include ../config.mk

.PHONY: default
default: $(addsuffix .o,$(basename $(wildcard *.cpp *.F)))

%.o: %.cpp
	$(CXX) -c $(CFLAGS) $<

%.o: %.F
	$(FC) -c $(FFLAGS) $<

clean:
	@rm -f *.o *.mod
