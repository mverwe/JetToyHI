
# Makefile generated automatically by scripts/mkcxx.pl '-a' 'test' '-f' '-s' '-1' '-r' '-8'
# run 'make make' to update it if you add new files

CXX = g++  # for macs - otherwise get c++ = clang
CXXFLAGS = -Wall -g -O2

# also arrange for fortran support
FC = gfortran
FFLAGS = -Wall -O2
CXXFLAGS += -std=c++11
LDFLAGS += -std=c++11

FJCONFIG = /afs/cern.ch/user/m/mverweij/work/soft/toy/fastjet-install/bin/fastjet-config
INCLUDE += `$(FJCONFIG) --cxxflags`
LIBRARIES  += `$(FJCONFIG) --libs --plugins` -lfastjetcontribfragile

PYTHIA8LOCATION = /afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67c/pythia8/226/x86_64-slc6-gcc48-opt
INCLUDE += -I$(PYTHIA8LOCATION)/include
LIBRARIES  += -L$(PYTHIA8LOCATION)/lib -lpythia8
LIBRARIES += -lgsl -lgslcblas -lm

INCLUDE += -I/usr/include


INCLUDE += `root-config --cflags`
LIBRARIES  += `root-config --libs`
INCLUDE += $(LCLINCLUDE)

COMMONSRC = 
F77SRC = 
COMMONOBJ = 

PROGSRC = runCreatePythiaEvents.cc runCreateThermalEvents.cc runtest.cc
PROGOBJ = runCreatePythiaEvents.o runCreateThermalEvents.o runtest.o

INCLUDE += 
LIBRARIES += 


all:  runCreatePythiaEvents runCreateThermalEvents runtest libtest.a


runCreatePythiaEvents: runCreatePythiaEvents.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runCreateThermalEvents: runCreateThermalEvents.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runtest: runtest.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

libtest.a: $(COMMONOBJ)
	ar cru libtest.a $(COMMONOBJ)
	ranlib libtest.a


make:
	scripts/mkcxx.pl '-a' 'test' '-f' '-s' '-1' '-r' '-8'

clean:
	rm -vf $(COMMONOBJ) $(PROGOBJ)

realclean: clean
	rm -vf  runCreatePythiaEvents runCreateThermalEvents runtest libtest.a

.cc.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.cpp.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.C.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.f.o:         $<
	$(FC) $(FFLAGS) -c $< -o $@
.f90.o:         $<
	$(FC) $(FFLAGS) -c $< -o $@


depend:
	makedepend  $(LCLINCLUDE) -Y --   -- $(COMMONSRC) $(PROGSRC)
# DO NOT DELETE

runCreatePythiaEvents.o: include/ProgressBar.h include/pythiaEvent.hh
runCreatePythiaEvents.o: include/extraInfo.hh
runCreateThermalEvents.o: include/ProgressBar.h include/thermalEvent.hh
runtest.o: include/ProgressBar.h include/jetCollection.hh
runtest.o: include/thermalEvent.hh include/pythiaEvent.hh
runtest.o: include/extraInfo.hh include/csSubtractor.hh
runtest.o: include/skSubtractor.hh include/softDropGroomer.hh
runtest.o: include/jetCollection.hh include/softDropCounter.hh
runtest.o: include/treeWriter.hh include/jetMatcher.hh include/randomCones.hh
