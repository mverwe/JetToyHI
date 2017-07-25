
# Makefile generated automatically by scripts/mkcxx.pl '-f' '-s' '-1' '-r' '-8' '-IPU14' '-l' '-LPU14 -lPU14 -lz'
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

PROGSRC = runCreatePythiaEvents.cc runCreateThermalEvents.cc runFromFile.cc runtest.cc
PROGOBJ = runCreatePythiaEvents.o runCreateThermalEvents.o runFromFile.o runtest.o

INCLUDE += 
LIBRARIES += -LPU14 -lPU14 -lz


all:  runCreatePythiaEvents runCreateThermalEvents runFromFile runtest 


runCreatePythiaEvents: runCreatePythiaEvents.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runCreateThermalEvents: runCreateThermalEvents.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runFromFile: runFromFile.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runtest: runtest.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)


make:
	scripts/mkcxx.pl '-f' '-s' '-1' '-r' '-8' '-IPU14' '-l' '-LPU14 -lPU14 -lz'

clean:
	rm -vf $(COMMONOBJ) $(PROGOBJ)

realclean: clean
	rm -vf  runCreatePythiaEvents runCreateThermalEvents runFromFile runtest 

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
runCreatePythiaEvents.o: include/extraInfo.hh include/extraInfo.hh
runCreatePythiaEvents.o: PU14/CmdLine.hh
runCreateThermalEvents.o: include/ProgressBar.h include/thermalEvent.hh
runCreateThermalEvents.o: PU14/CmdLine.hh
runFromFile.o: include/ProgressBar.h PU14/EventMixer.hh PU14/CmdLine.hh
runFromFile.o: PU14/EventSource.hh PU14/GSLRandom.hh PU14/CmdLine.hh
runFromFile.o: PU14/PU14.hh PU14/HepPID/ParticleIDMethods.hh
runFromFile.o: include/jetCollection.hh include/thermalEvent.hh
runFromFile.o: include/pythiaEvent.hh include/extraInfo.hh
runFromFile.o: include/csSubtractor.hh include/csSubtractorFullEvent.hh
runFromFile.o: include/skSubtractor.hh include/softDropGroomer.hh
runFromFile.o: include/jetCollection.hh include/treeWriter.hh
runFromFile.o: include/jetMatcher.hh include/randomCones.hh
runtest.o: PU14/CmdLine.hh include/ProgressBar.h include/jetCollection.hh
runtest.o: include/thermalEvent.hh include/pythiaEvent.hh
runtest.o: include/extraInfo.hh include/csSubtractor.hh
runtest.o: include/csSubtractorFullEvent.hh include/skSubtractor.hh
runtest.o: include/softDropGroomer.hh include/jetCollection.hh
runtest.o: include/softDropCounter.hh include/treeWriter.hh
runtest.o: include/jetMatcher.hh include/randomCones.hh
