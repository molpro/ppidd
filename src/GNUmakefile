
.PHONY: default
default: $(addsuffix .o,$(basename $(wildcard *.c *.F)))

%.o: %.c
	$(CC) -c $(CFLAGS) $<

%.o: %.F
	$(FC) -c $(FFLAGS) $<

clean:
	@rm -f *.o *.mod
