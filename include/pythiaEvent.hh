#ifndef pythiaEvent_h
#define pythiaEvent_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "Pythia8/Pythia.h"

//using namespace std;

//---------------------------------------------------------------
// Description
// This class generates a pythia8 event
// Author: M. Verweij
//---------------------------------------------------------------

class pythiaEvent {

private :
  Pythia8::Pythia pythia;
  double pthat_;
  unsigned int tune_;
  double rapMin_;
  double rapMax_;

public :
  pythiaEvent(double pthat = 120., unsigned int tune = 14, double rapMin = -3., double rapMax = 3.) :
    pthat_(pthat), tune_(tune), rapMin_(rapMin), rapMax_(rapMax)
  {
    
    // Generator. LHC process and output selection. Initialization.
    // tunes: http://home.thep.lu.se/~torbjorn/pythia82html/Tunes.html
    pythia.readString("Beams:eCM = 5002.");
    pythia.readString("HardQCD:all = on");
    pythia.readString(Form("PhaseSpace:pTHatMin = %.1f",pthat_));
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString(Form("Tune:pp = %d",tune_));
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 0");
    pythia.init();

  }

  std::vector<fastjet::PseudoJet> createPythiaEvent() {

    pythia.next();
    
    std::vector<fastjet::PseudoJet> particles;

    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].isFinal()) {
        fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
        if(p.rap()>rapMin_ && p.rap()<rapMax_)
          particles.push_back(p);
      }
    }

    return particles;
  }
};

#endif
