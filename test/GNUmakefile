
.PHONY: default
default: ppidd_test.exe ppidd_sf_test.exe ppidd_ctest.exe

ppidd_ctest.exe: ppidd_ctest.o
ifeq ($(BUILD),MPI2)
	$(CXX) $(CFLAGS) -o $@ $^ -L../lib -lppidd $(LIBS)
else
	$(FC) $(patsubst -Vaxlib,-nofor_main,$(FFLAGS)) -o $@ $^ -L../lib -lppidd $(LIBS)
endif

ppidd_test.exe: ppidd_test.o sizeofctypes.o sizeoffortypes.o ppidd_flush6.o
	$(FC) $(FFLAGS) -o $@ $^ -L../lib -lppidd $(LIBS)
ppidd_sf_test.exe: ppidd_sf_test.o ppidd_flush6.o
	$(FC) $(FFLAGS) -o $@ $^ -L../lib -lppidd $(LIBS)

ppidd_test.o: ppidd_test.F
ppidd_sf_test.o: ppidd_sf_test.F

%.o: %.cpp
	$(CXX) -c -I ../src $(CFLAGS) $<

%.o: %.F
	$(FC) -c -I ../src $(FFLAGS) $<

clean:
	@rm -f *.o *.exe
