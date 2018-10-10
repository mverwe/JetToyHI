
# Makefile generated automatically by ../scripts/mkcxx.pl '-a' 'PU14' '-f' '-s' '-1'
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
LIBRARIES += 
INCLUDE += 
INCLUDE += $(LCLINCLUDE)

COMMONSRC = CmdLine.cc EventMixer.cc EventSource.cc helpers.cc ParticleIDMethods.cc ParticleName.cc PU14.cc Version.cc zfstream.cc
F77SRC = 
COMMONOBJ = CmdLine.o EventMixer.o EventSource.o helpers.o ParticleIDMethods.o ParticleName.o PU14.o Version.o zfstream.o

PROGSRC = 
PROGOBJ = 

INCLUDE += 
LIBRARIES += 


all:  libPU14.a


libPU14.a: $(COMMONOBJ)
	ar cru libPU14.a $(COMMONOBJ)
	ranlib libPU14.a


make:
	../scripts/mkcxx.pl '-a' 'PU14' '-f' '-s' '-1'

clean:
	rm -vf $(COMMONOBJ) $(PROGOBJ)

realclean: clean
	rm -vf  libPU14.a

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
