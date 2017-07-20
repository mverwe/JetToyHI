#ifndef pythiaEvent_h
#define pythiaEvent_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "Pythia8/Pythia.h"

//using namespace std;

class pythiaEvent {

private :
  Pythia8::Pythia pythia;

public :
  pythiaEvent() {

    // Generator. LHC process and output selection. Initialization.
    // MV: selected some ATLAS tune - no idea how good it is - but missing LHPDF so running default tune now
    pythia.readString("Beams:eCM = 5002.");
    pythia.readString("HardQCD:all = on");
    pythia.readString("PhaseSpace:pTHatMin = 200.");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    //pythia.readString("Tune:pp = 19");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 0");
    pythia.init();

  }

  std::vector<fastjet::PseudoJet> createPythiaEvent() {

    pythia.next();
    
    std::vector<fastjet::PseudoJet> particles;

    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].isFinal()) {
        particles.push_back(fastjet::PseudoJet(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e()));
      }
    }

    return particles;
  }
};

#endif
