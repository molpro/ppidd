include ../config.mk

ifeq ($(FC),)
SOURCEFILES+=$(filter-out %fortran.cpp,$(wildcard *.cpp))
else
SOURCEFILES+=$(wildcard *.cpp *.F)
endif

.PHONY: default
default: $(addsuffix .o,$(basename $(SOURCEFILES)))

%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $<

%.o: %.F
	$(FC) -c $(FFLAGS) $<

clean:
	@rm -f *.o *.mod
