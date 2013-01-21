
.PHONY: default
default: ppidd_test.exe ppidd_sf_test.exe ppidd_ctest.exe

ppidd_ctest.exe: ppidd_ctest.o
ifeq ($(BUILD),MPI2)
	$(CC) $(CFLAGS) -o $@ $^ -L../lib -lppidd $(LIBS)
else
	$(FC) $(patsubst -Vaxlib,-nofor_main,$(FFLAGS)) -o $@ $^ -L../lib -lppidd $(LIBS)
endif

ppidd_test.exe: ppidd_test.o sizeofctypes.o sizeoffortypes.o ppidd_init_mod.o ppidd_flush6.o
	$(FC) $(FFLAGS) -o $@ $^ -L../lib -lppidd $(LIBS)
ppidd_sf_test.exe: ppidd_sf_test.o ppidd_init_mod.o ppidd_flush6.o
	$(FC) $(FFLAGS) -o $@ $^ -L../lib -lppidd $(LIBS)

ppidd_test.o: ppidd_test.F ppidd_init_mod.mod
ppidd_sf_test.o: ppidd_sf_test.F ppidd_init_mod.mod

%.o: %.c
	$(CC) -c -I ../src $(CFLAGS) $<

%.o %.mod: %.F
	$(FC) -c $(FFLAGS) $<

clean:
	@rm -f *.mod *.o *.exe
