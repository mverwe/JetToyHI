#ifndef thermalAlice_h
#define thermalAlice_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <random>

#include "extraInfo.hh"

//ROOT stuff
#include <TRandom3.h>
#include <TMath.h>
#include "TF1.h"

using namespace fastjet;
using namespace std;

//---------------------------------------------------------------
// Description
// This class generates a thermal event following Gamma distribution
// Based on https://github.com/ezradlesser/pyjetty/blob/28edacca74f7be4d7c22afd820d789ecfe8abe6d/pyjetty/alice_analysis/process/base/thermal_generator.py
//---------------------------------------------------------------

class thermalAlice {

private :
  TF1              *funcThrm_;

public :  
  std::vector<fastjet::PseudoJet> createThermalEventAlice() {

    std::vector<fastjet::PseudoJet> particles;

    double meanN = 2500.;
    double sigmaN = 500.;
    std::random_device rd {};
    std::mt19937 gen {rd()};
    std::normal_distribution<> d {meanN, sigmaN}; // mean , sigma
    
    unsigned int Nparticles = std::round(d(gen));

    for(unsigned int i = 0; i<Nparticles; ++i) {
      std::random_device rd_pt;
      std::mt19937 gen_pt(rd_pt());
      std::gamma_distribution<> d_pt(2, 0.3675);

      //pt from gamam function
      double pt = d_pt(gen_pt);
      //random phi

      double phi = gRandom->Rndm() * TMath::TwoPi();;
      //random rapidity

      double rap = gRandom->Rndm() * 1.8 - 0.9; // only positiveee need also negative
      //pion mass

      double mass = 0.1395; // pion mass
      
      fastjet::PseudoJet p4;
      p4.reset_momentum_PtYPhiM(pt,rap,phi,mass);
      particles.push_back(p4);
    }
    return particles;
  }

};

#endif
