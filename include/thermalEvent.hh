#ifndef thermalEvent_h
#define thermalEvent_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "extraInfo.hh"

//ROOT stuff
#include <TRandom3.h>
#include <TMath.h>
#include "TF1.h"

using namespace fastjet;
using namespace std;

//---------------------------------------------------------------
// Description
// This class generates a thermal event following Boltzman distribution
// Author: M. Verweij
//---------------------------------------------------------------

class thermalEvent {

private :
  double            meanpt_;
  unsigned int      mult_;
  double            rapMin_;
  double            rapMax_;
  double            neutralFrac_;
  TF1              *funcThrm_;

public :
  thermalEvent(unsigned int mult = 12000, double meanpt = 0.7, double rapMin = -3., double rapMax = 3.,double neutralFrac = 0.33) :
    meanpt_(meanpt),
    mult_(mult),
    rapMin_(rapMin),
    rapMax_(rapMax),
    neutralFrac_(neutralFrac)
  {
    
    funcThrm_ = new TF1("funcThrm_","[0]*TMath::Power([1], 2)*x*TMath::Exp(-[1]*x)", 0.2, 200.);
    funcThrm_->SetParNames("Amplitude", "b (GeV/c)^{-1}");
    funcThrm_->SetParameters(1.,2./meanpt_);
  }

  void setMult(unsigned int m) { mult_ = m; }
  void setMult(double mpt)     { meanpt_ = mpt; }
  void setRapidityRange(double min, double max) { rapMin_ =  min; rapMax_ =  max; }
  
  std::vector<fastjet::PseudoJet> createThermalEvent() {

    std::vector<fastjet::PseudoJet> particles;

    for(unsigned int i = 0; i<mult_; ++i) {
      double pt = funcThrm_->GetRandom();
      //random phi
      double phimin = 0.;
      double phimax = TMath::TwoPi();
      double phi = gRandom->Rndm() * (phimax - phimin) + phimin;
      //random rapidity
      double rap = gRandom->Rndm() * (rapMax_ - rapMin_) + rapMin_;

      int pdgid = 22;
      double mass = 0.0;
      if(gRandom->Rndm()>neutralFrac_) {
        mass = 0.1395;
        int charge = 1;
        if(gRandom->Rndm()<0.5) charge = -1;
        pdgid = charge*211;
      }
      //std::cout << "mass: " << mass << std::endl;
      
      fastjet::PseudoJet p4;
      p4.reset_momentum_PtYPhiM(pt,rap,phi,mass);
      p4.set_user_info(new extraInfo(pdgid, 1));
      
      particles.push_back(p4);
    }
    return particles;
  }
};

#endif
