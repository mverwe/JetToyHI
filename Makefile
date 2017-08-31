
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

PROGSRC = runCreatePythiaEvents.cc runCreatePythiaEventsPartonLevel.cc runCreateThermalEvents.cc runCSVariations.cc runFromFile.cc runJetPerformance.cc runJewelSub.cc runSDGenVarious.cc runtest.cc
PROGOBJ = runCreatePythiaEvents.o runCreatePythiaEventsPartonLevel.o runCreateThermalEvents.o runCSVariations.o runFromFile.o runJetPerformance.o runJewelSub.o runSDGenVarious.o runtest.o

INCLUDE += 
LIBRARIES += -LPU14 -lPU14 -lz


all:  runCreatePythiaEvents runCreatePythiaEventsPartonLevel runCreateThermalEvents runCSVariations runFromFile runJetPerformance runJewelSub runSDGenVarious runtest 


runCreatePythiaEvents: runCreatePythiaEvents.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runCreatePythiaEventsPartonLevel: runCreatePythiaEventsPartonLevel.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runCreateThermalEvents: runCreateThermalEvents.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runCSVariations: runCSVariations.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runFromFile: runFromFile.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runJetPerformance: runJetPerformance.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runJewelSub: runJewelSub.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runSDGenVarious: runSDGenVarious.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runtest: runtest.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)


make:
	scripts/mkcxx.pl '-f' '-s' '-1' '-r' '-8' '-IPU14' '-l' '-LPU14 -lPU14 -lz'

clean:
	rm -vf $(COMMONOBJ) $(PROGOBJ)

realclean: clean
	rm -vf  runCreatePythiaEvents runCreatePythiaEventsPartonLevel runCreateThermalEvents runCSVariations runFromFile runJetPerformance runJewelSub runSDGenVarious runtest 

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
runCreatePythiaEventsPartonLevel.o: include/ProgressBar.h
runCreatePythiaEventsPartonLevel.o: include/pythiaEvent.hh
runCreatePythiaEventsPartonLevel.o: include/extraInfo.hh include/extraInfo.hh
runCreatePythiaEventsPartonLevel.o: PU14/CmdLine.hh
runCreateThermalEvents.o: include/ProgressBar.h include/thermalEvent.hh
runCreateThermalEvents.o: include/extraInfo.hh PU14/CmdLine.hh
runCSVariations.o: include/ProgressBar.h PU14/EventMixer.hh PU14/CmdLine.hh
runCSVariations.o: PU14/EventSource.hh PU14/CmdLine.hh PU14/PU14.hh
runCSVariations.o: PU14/HepPID/ParticleIDMethods.hh include/jetCollection.hh
runCSVariations.o: include/csSubtractor.hh PU14/PU14.hh
runCSVariations.o: include/csSubtractorFullEvent.hh include/skSubtractor.hh
runCSVariations.o: include/softDropGroomer.hh include/jetCollection.hh
runCSVariations.o: include/treeWriter.hh include/jetMatcher.hh
runCSVariations.o: include/randomCones.hh include/Angularity.hh
runFromFile.o: include/ProgressBar.h PU14/EventMixer.hh PU14/CmdLine.hh
runFromFile.o: PU14/EventSource.hh PU14/CmdLine.hh PU14/PU14.hh
runFromFile.o: PU14/HepPID/ParticleIDMethods.hh include/jetCollection.hh
runFromFile.o: include/csSubtractor.hh PU14/PU14.hh
runFromFile.o: include/csSubtractorFullEvent.hh include/skSubtractor.hh
runFromFile.o: include/softDropGroomer.hh include/jetCollection.hh
runFromFile.o: include/treeWriter.hh include/jetMatcher.hh
runFromFile.o: include/randomCones.hh include/Angularity.hh
runFromFile.o: include/jewelMatcher.hh
runJetPerformance.o: include/ProgressBar.h PU14/EventMixer.hh PU14/CmdLine.hh
runJetPerformance.o: PU14/EventSource.hh PU14/CmdLine.hh PU14/PU14.hh
runJetPerformance.o: PU14/HepPID/ParticleIDMethods.hh
runJetPerformance.o: include/jetCollection.hh include/csSubtractor.hh
runJetPerformance.o: PU14/PU14.hh include/csSubtractorFullEvent.hh
runJetPerformance.o: include/skSubtractor.hh include/softDropGroomer.hh
runJetPerformance.o: include/jetCollection.hh include/treeWriter.hh
runJetPerformance.o: include/jetMatcher.hh include/randomCones.hh
runJetPerformance.o: include/Angularity.hh
runJewelSub.o: include/ProgressBar.h PU14/EventMixer.hh PU14/CmdLine.hh
runJewelSub.o: PU14/EventSource.hh PU14/CmdLine.hh PU14/PU14.hh
runJewelSub.o: PU14/HepPID/ParticleIDMethods.hh include/jetCollection.hh
runJewelSub.o: include/csSubtractor.hh PU14/PU14.hh
runJewelSub.o: include/csSubtractorFullEvent.hh include/skSubtractor.hh
runJewelSub.o: include/softDropGroomer.hh include/jetCollection.hh
runJewelSub.o: include/treeWriter.hh include/jetMatcher.hh
runJewelSub.o: include/randomCones.hh include/Angularity.hh
runJewelSub.o: include/jewelMatcher.hh
runSDGenVarious.o: include/ProgressBar.h PU14/EventMixer.hh PU14/CmdLine.hh
runSDGenVarious.o: PU14/EventSource.hh PU14/CmdLine.hh PU14/PU14.hh
runSDGenVarious.o: PU14/HepPID/ParticleIDMethods.hh include/jetCollection.hh
runSDGenVarious.o: include/csSubtractor.hh PU14/PU14.hh
runSDGenVarious.o: include/csSubtractorFullEvent.hh include/skSubtractor.hh
runSDGenVarious.o: include/softDropGroomer.hh include/jetCollection.hh
runSDGenVarious.o: include/treeWriter.hh include/jetMatcher.hh
runSDGenVarious.o: include/randomCones.hh include/Angularity.hh
runtest.o: PU14/CmdLine.hh include/ProgressBar.h include/jetCollection.hh
runtest.o: include/thermalEvent.hh include/extraInfo.hh
runtest.o: include/pythiaEvent.hh include/csSubtractor.hh PU14/PU14.hh
runtest.o: include/csSubtractorFullEvent.hh include/skSubtractor.hh
runtest.o: include/softDropGroomer.hh include/jetCollection.hh
runtest.o: include/softDropCounter.hh include/treeWriter.hh
runtest.o: include/jetMatcher.hh include/randomCones.hh
