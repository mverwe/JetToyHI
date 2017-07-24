#ifndef randomCones_h
#define randomCones_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "fastjet/contrib/ConstituentSubtractor.hh"

#include "TVector2.h"
#include "TRandom3.h"
#include "TMath.h"

using namespace std;

//---------------------------------------------------------------
// Description
// This class runs random cones on a event
// and stores the random cones without subtraction in a PseudoJet vector
// Author: M. Verweij
//---------------------------------------------------------------

class randomCones {

private :
  unsigned int nCones_;
  double rParam_;
  double etaMax_;
  std::vector<fastjet::PseudoJet> fjInputs_;
  TRandom3 *rnd_;

  double deltaR(const double phi1, const double phi2, const double eta1, const double eta2);

public :
  randomCones(unsigned int nCones = 4, double rParam = 0.4, double etaMax = 2.3);

  void setInputParticles(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }
  
  std::vector<fastjet::PseudoJet> run();
};

randomCones::randomCones(unsigned int nCones, double rParam, double etaMax)
  : nCones_(nCones),
    rParam_(rParam),
    etaMax_(etaMax),
    rnd_()
{
  if(!rnd_) rnd_ = new TRandom3(0);
}


std::vector<fastjet::PseudoJet>  randomCones::run() {
  
  if(!rnd_) rnd_ = new TRandom3(0);
    
  std::vector<fastjet::PseudoJet> cones;
  cones.reserve(nCones_);

  double minPhi = 0.;
  double maxPhi = 2.*TMath::Pi();
    
  for(unsigned int i = 0; i<nCones_; ++i) {

    //pick random position for random cone
    double etaRC = rnd_->Rndm() * (etaMax_ - -1.*etaMax_) + -1.*etaMax_;
    double phiRC = rnd_->Rndm() * (maxPhi - minPhi) + minPhi;

    double ptRC = 0.;
    for(fastjet::PseudoJet part : fjInputs_) {
      double dr = deltaR(part.phi(),phiRC,part.eta(),etaRC);
      if(dr<rParam_) {
        //std::cout << "phi: " << part.phi() << " " << phiRC << " eta: " << part.eta() << " " << etaRC << std::endl;
        ptRC+=part.pt();
      }
    }

    fastjet::PseudoJet cone;
    cone.reset_momentum_PtYPhiM(ptRC,etaRC,phiRC,0.);
    cones.push_back(cone);
  }
    
  return cones;
}

double randomCones::deltaR(const double phi1, const double phi2, const double eta1, const double eta2) {
  //calculate distance
  double dPhi = phi1 - phi2;
  double dEta = eta1 - eta2;
  dPhi = TVector2::Phi_mpi_pi(dPhi);
  return TMath::Sqrt(dPhi * dPhi + dEta * dEta);
}

#endif
