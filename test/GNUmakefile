include ../config.mk

.PHONY: default
ifeq ($(FC),)
default: ppidd_ctest.exe
else
default: ppidd_test.exe ppidd_sf_test.exe ppidd_ctest.exe
endif

ppidd_ctest.exe: ppidd_ctest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ -L.. -lppidd $(LIBS)

ppidd_test.exe: ppidd_test.o sizeofctypes.o sizeoffortypes.o ppidd_flush6.o
	$(FC) $(FFLAGS) -o $@ $^ -L.. -lppidd $(LIBS)

ppidd_sf_test.exe: ppidd_sf_test.o ppidd_flush6.o
	$(FC) $(FFLAGS) -o $@ $^ -L.. -lppidd $(LIBS)

ppidd_test.o: ppidd_test.F
ppidd_sf_test.o: ppidd_sf_test.F

%.o: %.cpp
	$(CXX) $(CPPFLAGS) -c -I ../src $(CXXFLAGS) $<

%.o: %.F
	$(FC) $(CPPFLAGS) -c -I ../src $(FCFLAGS) $<
