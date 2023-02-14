
# Makefile generated automatically by ../scripts/mkcxx.pl '-a' 'PU14' '-f' '-s' '-1'
# run 'make make' to update it if you add new files

CXX = g++  # for macs - otherwise get c++ = clang
CXXFLAGS = -Wall -g -O2

# also arrange for fortran support
FC = gfortran
FFLAGS = -Wall -O2
CXXFLAGS += -std=c++11
LDFLAGS += -std=c++11

FJCONFIG = /Users/mverweij/soft/fastjet-3.4.0/../fastjet340-install/bin/fastjet-config
FJLOC = /Users/mverweij/soft/fastjet-3.4.0/../fastjet340-install

INCLUDE += `$(FJCONFIG) --cxxflags`
LIBRARIES  += -L$(FJLOC)/lib -lRecursiveTools `$(FJCONFIG) --libs --plugins` -lfastjetcontribfragile
LIBRARIES += -L/usr/local/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas

INCLUDE += -I/usr/local/Cellar/gsl/2.7.1/include

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
# DO NOT DELETE

CmdLine.o: CmdLine.hh
EventMixer.o: EventMixer.hh CmdLine.hh EventSource.hh PU14.hh
EventMixer.o: HepPID/ParticleIDMethods.hh helpers.hh
EventSource.o: EventSource.hh FastIStringStream.hh PU14.hh
EventSource.o: HepPID/ParticleIDMethods.hh zfstream.h
helpers.o: helpers.hh
ParticleIDMethods.o: HepPID/ParticleIDMethods.hh HepPID/ParticleName.hh
ParticleName.o: HepPID/ParticleName.hh HepPID/ParticleIDMethods.hh
ParticleName.o: HepPID/Version.hh
PU14.o: PU14.hh HepPID/ParticleIDMethods.hh
Version.o: HepPID/Version.hh
zfstream.o: zfstream.h
