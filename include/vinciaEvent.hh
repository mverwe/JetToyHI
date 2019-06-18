#ifndef pythiaEvent_h
#define pythiaEvent_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "Pythia8/Pythia.h"
#include "Vincia/Vincia.h"

#include "extraInfo.hh"

//using namespace std;
using namespace Pythia8;
using namespace Vincia;

//---------------------------------------------------------------
// Description
// This class generates a vincia event
// Author: M. Verweij
//---------------------------------------------------------------

class vinciaEvent {

private :
  Pythia8::Pythia pythia;
  //VinciaPlugin vincia(&pythia);
  //, "/Users/mverweij/soft/pythia8240/vincia-2.2.04/share/Vincia/xmldoc/VinciaIndex.xml");
  //
  double pthat_;
  unsigned int tune_;
  double rapMin_;
  double rapMax_;
  bool   partonLevel_;

  std::vector<fastjet::PseudoJet> partons;

public :
  vinciaEvent(double pthat = 120., unsigned int tune = 14, double rapMin = -3., double rapMax = 3., bool partonLevel = false, std::string xmlDocLoc = "/Users/mverweij/soft/pythia8240/vincia-2.2.04/share/Vincia/xmldoc");
  std::vector<fastjet::PseudoJet> createVinciaEvent();
  
  std::vector<fastjet::PseudoJet> getPartonList() const { return partons; }

};

vinciaEvent::vinciaEvent(double pthat, unsigned int tune, double rapMin, double rapMax, bool partonLevel, std::string xmlDocLoc) :
  pthat_(pthat), tune_(tune), rapMin_(rapMin), rapMax_(rapMax), partonLevel_(partonLevel)
{
  VinciaPlugin vincia(&pythia, xmlDocLoc);
  
  // Generator. LHC process and output selection. Initialization.
  // tunes: http://home.thep.lu.se/~torbjorn/pythia82html/Tunes.html
  pythia.readString("Beams:eCM = 5002.");
  pythia.readString("HardQCD:all = on");
  //pythia.readString("HardQCD:gg2gg = on");
  pythia.readString("Vincia:FSR = on");
  pythia.readString("Vincia:ISR = on");
  //pythia.readString("Vincia:biasAll = 5");
  //pythia.readString("Vincia:biasBottom = 50");

  //pythia.readString("Vincia:nGluonToQuarkF = 5");
  // pythia.readString("Vincia:QQSplitFF:chargeFactor = 1.0");
  
  pythia.readString(Form("PhaseSpace:pTHatMin = %.1f",pthat_));
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Switch on/off MPI and hadronization
  pythia.readString("PartonLevel:MPI = on");
     
  pythia.readString(Form("Tune:pp = %d",tune_));
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  if(partonLevel_) {
    pythia.readString("HadronLevel:all = off");
  }
    
  vincia.init();
  //pythia.init();
    
  Printf("\n");

}

std::vector<fastjet::PseudoJet> vinciaEvent::createVinciaEvent() {

  pythia.next();
    
  std::vector<fastjet::PseudoJet> particles;
  partons.clear();

  // Printf("\n\n new event");

  for (int i = 0; i < pythia.event.size(); ++i) {
    if (pythia.event[i].isFinal()) {
      fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
      p.set_user_info(new extraInfo(pythia.event[i].id(), 0)); 
      if(p.rap()>rapMin_ && p.rap()<rapMax_)
        particles.push_back(p);
    } else if (pythia.event[i].status()==-23) {
      //Printf("Initial parton: %d, pT: %f  eta: %f  phi: %f  mass: %f",pythia.event[i].id(),pythia.event[i].pT(),pythia.event[i].eta(),pythia.event[i].phi(),pythia.event[i].m());
      // if(pythia.event[i].daughter1()>0)
      //   Printf("Daughter 1: %d, pT: %f  eta: %f  phi: %f",pythia.event[pythia.event[i].daughter1()].id(),pythia.event[pythia.event[i].daughter1()].pT(),pythia.event[pythia.event[i].daughter1()].eta(),pythia.event[pythia.event[i].daughter1()].phi());
      // if(pythia.event[i].daughter2()>0)
      //   Printf("Daughter 2: %d, pT: %f  eta: %f  phi: %f",pythia.event[pythia.event[i].daughter2()].id(),pythia.event[pythia.event[i].daughter2()].pT(),pythia.event[pythia.event[i].daughter2()].eta(),pythia.event[pythia.event[i].daughter2()].phi());

      fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
      p.set_user_info(new extraInfo(pythia.event[i].id(), -1)); 
      partons.push_back(p);
    }
    else if(pythia.event[i].status()==-51 && abs(pythia.event[i].id())==5) {
      //check for b quarks
      // Printf("!!!!!!!!!!!!!! %d status: %d Found b quark from parton branching: %d, pT: %f  eta: %f  phi: %f  mass: %f",i,pythia.event[i].status(),pythia.event[i].id(),pythia.event[i].pT(),pythia.event[i].eta(),pythia.event[i].phi(),pythia.event[i].m());
    }
    else if(abs(pythia.event[i].id())==5) {
      //check for a b quarks
      // Printf("!!!!!!!!!!!!!! %d status: %d Found a b quark: %d, pT: %f  eta: %f  phi: %f  mass: %f",i,pythia.event[i].status(),pythia.event[i].id(),pythia.event[i].pT(),pythia.event[i].eta(),pythia.event[i].phi(),pythia.event[i].m());
    }
    // else if(pythia.event[i].status()==-33) {// && abs(pythia.event[i].id())==5) {
    //   Printf("!!!!!!!!!!!!!! %d Found outgoing subsequent subprocess: %d, pT: %f  eta: %f  phi: %f",i,pythia.event[i].id(),pythia.event[i].pT(),pythia.event[i].eta(),pythia.event[i].phi());
    // }

    // if(pythia.event[i].status()<-20 && pythia.event[i].status()>-60) {
    //   Printf("!!!!!!!!!!!!!! %d Found particle with status between 20 and 60: %d, pT: %f  eta: %f  phi: %f",i,pythia.event[i].id(),pythia.event[i].pT(),pythia.event[i].eta(),pythia.event[i].phi());
    // }

    //Printf("!!!!!!!!!!!!!! %d Found particle with status between 20 and 60: %d, pT: %f  eta: %f  phi: %f",i,pythia.event[i].id(),pythia.event[i].pT(),pythia.event[i].eta(),pythia.event[i].phi());


  }

  // pythia.event.list();
    
  //  pythia.stat();
    
  return particles;
}

#endif
