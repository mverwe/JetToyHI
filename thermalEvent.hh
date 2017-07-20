#ifndef thermalEvent_h
#define thermalEvent_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

//ROOT stuff
#include <TRandom3.h>
#include <TMath.h>
#include "TF1.h"

using namespace fastjet;
using namespace std;

class thermalEvent {

private :
  double            meanpt_;
  unsigned int      mult_;
  double            etaMin_;
  double            etaMax_;
  TF1              *funcThrm_;

public :
  thermalEvent() {

    meanpt_ = 0.7;
    mult_   = 12000;
    etaMin_ = -3.;
    etaMax_ = 3.;
    
    funcThrm_ = new TF1("funcThrm_","[0]*TMath::Power([1], 2)*x*TMath::Exp(-[1]*x)", 0.2, 200.);
    funcThrm_->SetParNames("Amplitude", "b (GeV/c)^{-1}");
    funcThrm_->SetParameters(1.,2./meanpt_);
    
  }

  void setMult(unsigned int m) { mult_ = m; }
  void setMult(double mpt)     { meanpt_ = mpt; }
  void setEtaRange(double min, double max) { etaMin_ =  min; etaMax_ =  max; }
  
  std::vector<fastjet::PseudoJet> createThermalEvent() {

    std::vector<fastjet::PseudoJet> particles;

    for(unsigned int i = 0; i<mult_; ++i) {
      double pt = funcThrm_->GetRandom();
      //random phi
      double phimin = 0.;
      double phimax = TMath::TwoPi();
      double phi = gRandom->Rndm() * (phimax - phimin) + phimin;
      //random eta
      double eta = gRandom->Rndm() * (etaMax_ - etaMin_) + etaMin_;
      
      fastjet::PseudoJet p4;
      p4.reset_momentum_PtYPhiM(pt,eta,phi,0.1395);
      
      particles.push_back(p4);
    }
    return particles;
  }
};

#endif
