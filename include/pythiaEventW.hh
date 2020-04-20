#ifndef pythiaEventW_h
#define pythiaEventW_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "Pythia8/Pythia.h"

#include "extraInfo.hh"

//using namespace std;

//---------------------------------------------------------------
// Description
// This class generates a pythia8 event
// Author: M. Verweij
//---------------------------------------------------------------

class pythiaEventW {

private :
  Pythia8::Pythia pythia;
  double pthat_;
  unsigned int tune_;
  double rapMin_;
  double rapMax_;
  bool   partonLevel_;

  std::vector<fastjet::PseudoJet> partons;

public :
  pythiaEventW(double pthat = 120., unsigned int tune = 14, double rapMin = -3., double rapMax = 3., bool partonLevel = false);
  std::vector<fastjet::PseudoJet> createPythiaEvent();
  
  std::vector<fastjet::PseudoJet> getPartonList() const { return partons; }

};
  
pythiaEventW::pythiaEventW(double pthat, unsigned int tune, double rapMin, double rapMax, bool partonLevel) :
  pthat_(pthat), tune_(tune), rapMin_(rapMin), rapMax_(rapMax), partonLevel_(partonLevel)
{
    
  // Generator. LHC process and output selection. Initialization.
  // tunes: http://home.thep.lu.se/~torbjorn/pythia82html/Tunes.html
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 5002.");
  //pythia.readString("HardQCD:all = on");
  pythia.readString(Form("PhaseSpace:pTHatMin = %.1f",pthat_));
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString(Form("Tune:pp = %d",tune_));
  pythia.readString("Random:setSeed = on");
  pythia.readString("WeakDoubleBoson:ffbar2WW=on");
  pythia.readString("24:onIfAny = 1 2 3 4 5");
  pythia.readString("310:mayDecay = off");
  pythia.readString("3122:mayDecay = off");
  pythia.readString("3112:mayDecay = off");
  pythia.readString("3212:mayDecay = off");
  pythia.readString("3222:mayDecay = off");
  pythia.readString("3312:mayDecay = off");
  pythia.readString("3322:mayDecay = off");
  pythia.readString("3334:mayDecay = off");

  pythia.readString("Random:seed = 0");
  if(partonLevel_) {
    pythia.readString("HadronLevel:all = off");
  }
 
    
  pythia.init();
  }

std::vector<fastjet::PseudoJet> pythiaEventW::createPythiaEvent() {

  pythia.next(); //generate next event
    
  std::vector<fastjet::PseudoJet> particles;
  partons.clear(); //empty list before storing partons of new event

  pythia.event.list();
  
  for (int i = 0; i < pythia.event.size(); ++i) {
    if (pythia.event[i].isFinal()) {
      fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
      p.set_user_info(new extraInfo(pythia.event[i].id(), 0)); 
      if(p.rap()>rapMin_ && p.rap()<rapMax_)
        particles.push_back(p);
    } else if (pythia.event[i].status()==-23) {
      fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
      p.set_user_info(new extraInfo(pythia.event[i].id(), -1)); 
      partons.push_back(p);
    }

    //if(abs(pythia.event[i].id())<7) {
    //  std::cout << "px: " << pythia.event[i].px() << "  py: " << pythia.event[i].py() << "  pz: " << pythia.event[i].pz() << "  status: " << pythia.event[i].status() << "  id: " << pythia.event[i].id() << std::endl;
    //}

  }

  return particles;
}


#endif
