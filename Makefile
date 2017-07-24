SOURCES := $(wildcard *.cc)
EXECUTABLE := $(patsubst %.cc, %, $(SOURCES))

DEPS=include/thermalEvent.hh include/pythiaEvent.hh include/csSubtractor.hh include/skSubtractor.hh include/softDropGroomer.hh include/treeWriter.hh include/jetMatcher.hh include/randomCones.hh
# FASTJET=/Users/yichen/Programs/fastjet/install321-1026
# PYTHIA8LOCATION=/Users/yichen/Programs/pythia/pythia8219
FASTJET=/afs/cern.ch/user/m/mverweij/work/soft/toy/fastjet-install
PYTHIA8LOCATION=/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67c/pythia8/226/x86_64-slc6-gcc48-opt
LIBDIRARCH=lib

CC=g++
CFLAGS=-c -I$(ROOTSYS)/include #-I$(PYTHIA8LOCATION)/include 

# LDFLAGS=-L$(ROOTSYS)/lib -L$(PYTHIA8LOCATION)/lib -lpythia8 $(shell root-config --glibs) $(shell root-config --cflags) -I$(PYTHIA8LOCATION)/include `$(FASTJET)/bin/fastjet-config --cxxflags --plugins` $(FASTJET)/lib/*.a -L$(PYTHIA8LOCATION)/lib -lpythia8
LDFLAGS=-L$(ROOTSYS)/lib -L$(PYTHIA8LOCATION)/lib -lpythia8 $(shell root-config --glibs) $(shell root-config --cflags) -I$(PYTHIA8LOCATION)/include `$(FASTJET)/bin/fastjet-config --cxxflags --libs --plugins` -lfastjetcontribfragile -L$(PYTHIA8LOCATION)/lib -lpythia8

all: $(EXECUTABLE)
$(EXECUTABLE): %: %.cc $(DEPS)
	$(CC) $(LDFLAGS) $< -o $@

.PHONY: clean all

clean:
	rm -rf *.o *~ $(EXECUTABLE)


