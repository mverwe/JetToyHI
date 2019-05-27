
# Makefile generated automatically by scripts/mkcxx.pl '-f' '-s' '-1' '-r' '-8' '-IPU14' '-l' '-LPU14 -lPU14 -lz'
# run 'make make' to update it if you add new files

CXX = g++  # for macs - otherwise get c++ = clang
CXXFLAGS = -Wall -g -O2

# also arrange for fortran support
FC = gfortran
FFLAGS = -Wall -O2
CXXFLAGS += -std=c++11
LDFLAGS += -std=c++11

FJCONFIG = /Users/mverweij/soft/fastjet330-install/bin/fastjet-config
INCLUDE += `$(FJCONFIG) --cxxflags`
LIBRARIES  += `$(FJCONFIG) --libs --plugins` -lfastjetcontribfragile

PYTHIA8LOCATION = /afs/cern.ch/user/l/lcunquei/pythia8235
INCLUDE += -I$(PYTHIA8LOCATION)/include
LIBRARIES  += -L$(PYTHIA8LOCATION)/lib -lpythia8
LIBRARIES += -L/usr/local/Cellar/gsl/2.5/lib -lgsl -lgslcblas

INCLUDE += -I/usr/local/Cellar/gsl/2.5/include


INCLUDE += `root-config --cflags`
LIBRARIES  += `root-config --glibs`
INCLUDE += $(LCLINCLUDE)

COMMONSRC = 
F77SRC = 
COMMONOBJ = 

PROGSRC = PU14toHepMC.cc runConversionQPYTHIA.cc runCreatePythiaEvents.cc runCreatePythiaEventsPartonLevel.cc runCreateThermalEvents.cc runCSVariations.cc runFromFile.cc runJetPerformance.cc runJetProfile.cc runJetProfileJewelSub.cc runJewelSub.cc runSDGenVarious.cc runSDGenVariousJewelSub.cc runSharedLayerSubtraction.cc runSimpleJetAnalysis.cc runtest.cc
PROGOBJ = PU14toHepMC.o runConversionQPYTHIA.o runCreatePythiaEvents.o runCreatePythiaEventsPartonLevel.o runCreateThermalEvents.o runCSVariations.o runFromFile.o runJetPerformance.o runJetProfile.o runJetProfileJewelSub.o runJewelSub.o runSDGenVarious.o runSDGenVariousJewelSub.o runSharedLayerSubtraction.o runSimpleJetAnalysis.o runtest.o

INCLUDE += 
LIBRARIES += -LPU14 -lPU14 -lz


all:  PU14toHepMC runConversionQPYTHIA runCreatePythiaEvents runCreatePythiaEventsPartonLevel runCreateThermalEvents runCSVariations runFromFile runJetPerformance runJetProfile runJetProfileJewelSub runJewelSub runSDGenVarious runSDGenVariousJewelSub runSharedLayerSubtraction runSimpleJetAnalysis runtest 


PU14toHepMC: PU14toHepMC.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runConversionQPYTHIA: runConversionQPYTHIA.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

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

runJetProfile: runJetProfile.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runJetProfileJewelSub: runJetProfileJewelSub.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runJewelSub: runJewelSub.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runSDGenVarious: runSDGenVarious.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runSDGenVariousJewelSub: runSDGenVariousJewelSub.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runSharedLayerSubtraction: runSharedLayerSubtraction.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runSimpleJetAnalysis: runSimpleJetAnalysis.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)

runtest: runtest.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)


make:
	scripts/mkcxx.pl '-f' '-s' '-1' '-r' '-8' '-IPU14' '-l' '-LPU14 -lPU14 -lz'

clean:
	rm -vf $(COMMONOBJ) $(PROGOBJ)

realclean: clean
	rm -vf  PU14toHepMC runConversionQPYTHIA runCreatePythiaEvents runCreatePythiaEventsPartonLevel runCreateThermalEvents runCSVariations runFromFile runJetPerformance runJetProfile runJetProfileJewelSub runJewelSub runSDGenVarious runSDGenVariousJewelSub runSharedLayerSubtraction runSimpleJetAnalysis runtest 

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

runConversionQPYTHIA.o: include/ProgressBar.h include/pythiaEvent.hh
runConversionQPYTHIA.o: include/extraInfo.hh include/extraInfo.hh
runConversionQPYTHIA.o: PU14/CmdLine.hh
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
runCSVariations.o: include/jewelMatcher.hh include/treeWriter.hh
runCSVariations.o: include/jetMatcher.hh include/randomCones.hh
runCSVariations.o: include/Angularity.hh
runFromFile.o: include/ProgressBar.h PU14/EventMixer.hh PU14/CmdLine.hh
runFromFile.o: PU14/EventSource.hh PU14/CmdLine.hh PU14/PU14.hh
runFromFile.o: PU14/HepPID/ParticleIDMethods.hh include/jetCollection.hh
runFromFile.o: include/csSubtractor.hh PU14/PU14.hh
runFromFile.o: include/csSubtractorFullEvent.hh include/skSubtractor.hh
runFromFile.o: include/softDropGroomer.hh include/jetCollection.hh
runFromFile.o: include/jewelMatcher.hh include/treeWriter.hh
runFromFile.o: include/jetMatcher.hh include/randomCones.hh
runFromFile.o: include/Angularity.hh include/jewelMatcher.hh
runJetPerformance.o: include/ProgressBar.h PU14/EventMixer.hh PU14/CmdLine.hh
runJetPerformance.o: PU14/EventSource.hh PU14/CmdLine.hh PU14/PU14.hh
runJetPerformance.o: PU14/HepPID/ParticleIDMethods.hh
runJetPerformance.o: include/jetCollection.hh include/csSubtractor.hh
runJetPerformance.o: PU14/PU14.hh include/csSubtractorFullEvent.hh
runJetPerformance.o: include/skSubtractor.hh include/softDropGroomer.hh
runJetPerformance.o: include/jetCollection.hh include/jewelMatcher.hh
runJetPerformance.o: include/treeWriter.hh include/jetMatcher.hh
runJetPerformance.o: include/randomCones.hh include/Angularity.hh
runJetProfile.o: include/ProgressBar.h PU14/EventMixer.hh PU14/CmdLine.hh
runJetProfile.o: PU14/EventSource.hh PU14/CmdLine.hh PU14/PU14.hh
runJetProfile.o: PU14/HepPID/ParticleIDMethods.hh include/jetCollection.hh
runJetProfile.o: include/csSubtractor.hh PU14/PU14.hh
runJetProfile.o: include/csSubtractorFullEvent.hh include/skSubtractor.hh
runJetProfile.o: include/softDropGroomer.hh include/jetCollection.hh
runJetProfile.o: include/jewelMatcher.hh include/treeWriter.hh
runJetProfile.o: include/jetMatcher.hh include/randomCones.hh
runJetProfile.o: include/Angularity.hh
runJetProfileJewelSub.o: include/ProgressBar.h PU14/EventMixer.hh
runJetProfileJewelSub.o: PU14/CmdLine.hh PU14/EventSource.hh PU14/CmdLine.hh
runJetProfileJewelSub.o: PU14/PU14.hh PU14/HepPID/ParticleIDMethods.hh
runJetProfileJewelSub.o: include/jetCollection.hh include/csSubtractor.hh
runJetProfileJewelSub.o: PU14/PU14.hh include/csSubtractorFullEvent.hh
runJetProfileJewelSub.o: include/skSubtractor.hh include/softDropGroomer.hh
runJetProfileJewelSub.o: include/jetCollection.hh include/jewelMatcher.hh
runJetProfileJewelSub.o: include/treeWriter.hh include/jetMatcher.hh
runJetProfileJewelSub.o: include/randomCones.hh include/Angularity.hh
runJewelSub.o: include/ProgressBar.h PU14/EventMixer.hh PU14/CmdLine.hh
runJewelSub.o: PU14/EventSource.hh PU14/CmdLine.hh PU14/PU14.hh
runJewelSub.o: PU14/HepPID/ParticleIDMethods.hh include/jetCollection.hh
runJewelSub.o: include/csSubtractor.hh PU14/PU14.hh
runJewelSub.o: include/csSubtractorFullEvent.hh include/skSubtractor.hh
runJewelSub.o: include/softDropGroomer.hh include/jetCollection.hh
runJewelSub.o: include/jewelMatcher.hh include/treeWriter.hh
runJewelSub.o: include/jetMatcher.hh include/randomCones.hh
runJewelSub.o: include/Angularity.hh include/jewelMatcher.hh
runSDGenVarious.o: include/ProgressBar.h PU14/EventMixer.hh PU14/CmdLine.hh
runSDGenVarious.o: PU14/EventSource.hh PU14/CmdLine.hh PU14/PU14.hh
runSDGenVarious.o: PU14/HepPID/ParticleIDMethods.hh include/jetCollection.hh
runSDGenVarious.o: include/csSubtractor.hh PU14/PU14.hh
runSDGenVarious.o: include/csSubtractorFullEvent.hh include/skSubtractor.hh
runSDGenVarious.o: include/softDropGroomer.hh include/jetCollection.hh
runSDGenVarious.o: include/jewelMatcher.hh include/treeWriter.hh
runSDGenVarious.o: include/jetMatcher.hh include/randomCones.hh
runSDGenVarious.o: include/Angularity.hh
runSDGenVariousJewelSub.o: include/ProgressBar.h PU14/EventMixer.hh
runSDGenVariousJewelSub.o: PU14/CmdLine.hh PU14/EventSource.hh
runSDGenVariousJewelSub.o: PU14/CmdLine.hh PU14/PU14.hh
runSDGenVariousJewelSub.o: PU14/HepPID/ParticleIDMethods.hh
runSDGenVariousJewelSub.o: include/jetCollection.hh include/csSubtractor.hh
runSDGenVariousJewelSub.o: PU14/PU14.hh include/csSubtractorFullEvent.hh
runSDGenVariousJewelSub.o: include/skSubtractor.hh include/softDropGroomer.hh
runSDGenVariousJewelSub.o: include/jetCollection.hh include/jewelMatcher.hh
runSDGenVariousJewelSub.o: include/treeWriter.hh include/jetMatcher.hh
runSDGenVariousJewelSub.o: include/randomCones.hh include/Angularity.hh
runSDGenVariousJewelSub.o: include/jewelMatcher.hh
runSharedLayerSubtraction.o: include/ProgressBar.h PU14/EventMixer.hh
runSharedLayerSubtraction.o: PU14/CmdLine.hh PU14/EventSource.hh
runSharedLayerSubtraction.o: PU14/CmdLine.hh PU14/PU14.hh
runSharedLayerSubtraction.o: PU14/HepPID/ParticleIDMethods.hh
runSharedLayerSubtraction.o: include/jetCollection.hh
runSharedLayerSubtraction.o: include/sharedLayerSubtractor.hh PU14/PU14.hh
runSharedLayerSubtraction.o: include/Angularity.hh include/treeWriter.hh
runSharedLayerSubtraction.o: include/jetCollection.hh include/jetMatcher.hh
runSimpleJetAnalysis.o: include/ProgressBar.h PU14/EventMixer.hh
runSimpleJetAnalysis.o: PU14/CmdLine.hh PU14/EventSource.hh PU14/CmdLine.hh
runSimpleJetAnalysis.o: PU14/PU14.hh PU14/HepPID/ParticleIDMethods.hh
runSimpleJetAnalysis.o: include/jetCollection.hh include/softDropGroomer.hh
runSimpleJetAnalysis.o: PU14/PU14.hh include/jetCollection.hh
runSimpleJetAnalysis.o: include/jewelMatcher.hh include/treeWriter.hh
runSimpleJetAnalysis.o: include/jetMatcher.hh include/Angularity.hh
runtest.o: PU14/CmdLine.hh include/ProgressBar.h include/jetCollection.hh
runtest.o: include/thermalEvent.hh include/extraInfo.hh
runtest.o: include/pythiaEvent.hh include/csSubtractor.hh PU14/PU14.hh
runtest.o: include/csSubtractorFullEvent.hh include/skSubtractor.hh
runtest.o: include/softDropGroomer.hh include/jetCollection.hh
runtest.o: include/jewelMatcher.hh include/softDropCounter.hh
runtest.o: include/treeWriter.hh include/jetMatcher.hh include/randomCones.hh
