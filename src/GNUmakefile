include ../config.mk

SOURCEFILES+=$(wildcard *.cpp)
ifneq ($(FC),)
SOURCEFILES+=$(wildcard *.F)
endif

.PHONY: default
default: $(addsuffix .o,$(basename $(SOURCEFILES)))

%.o: %.cpp
	$(CC) -c $(CFLAGS) $<

%.o: %.F
	$(FC) -c $(FFLAGS) $<

clean:
	@rm -f *.o *.mod
