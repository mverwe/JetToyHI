SOURCES=runtest.cc
DEPS=thermalEvent.hh pythiaEvent.hh csSubtractor.hh
EXECUTABLE=runtest
OBJECTS=$(SOURCES:.cpp=.o)
FASTJET=/afs/cern.ch/user/m/mverweij/work/soft/toy/fastjet-install
#FASTJET=/afs/cern.ch/sw/lcg/external/fastjet/3.2.0/x86_64-slc6-gcc48-opt/
#FASTJET=/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/fastjet/3.1.0
PYTHIA8LOCATION=/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67c/pythia8/226/x86_64-slc6-gcc48-opt
LIBDIRARCH=lib

CC=g++
CFLAGS=-c -I$(ROOTSYS)/include #-I$(PYTHIA8LOCATION)/include 
#LDFLAGS=$(shell root-config --glibs) $(shell root-config --cflags) `$(FASTJET)/bin/fastjet-config --cxxflags --libs --plugins`
#LDFLAGS=-L$(ROOTSYS)/lib $(shell root-config --glibs) $(shell root-config --cflags) `$(FASTJET)/bin/fastjet-config --cxxflags --libs --plugins` -L$(PYTHIA8LOCATION)/$(LIBDIRARCH) -lpythia8 -ldl

LDFLAGS=-L$(ROOTSYS)/lib -L$(PYTHIA8LOCATION)/lib -lpythia8 $(shell root-config --glibs) $(shell root-config --cflags) -I$(PYTHIA8LOCATION)/include `$(FASTJET)/bin/fastjet-config --cxxflags --libs --plugins` -lfastjetcontribfragile -L$(PYTHIA8LOCATION)/lib -lpythia8# `/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67c/pythia8/226/x86_64-slc6-gcc48-opt/bin/pythia8-config --cxxflags --libs` -lpythia8 #-L/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67c/pythia8/226/x86_64-slc6-gcc48-opt/lib -lpythia8

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(DEPS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o: $(PYTHIA8LOCATION)/libpythia8.a
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -rf *.o *~ $(EXECUTABLE)

