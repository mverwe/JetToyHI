#ifndef __jetProfile_HH__
#define __jetProfile_HH__

//------------------------------------------------------------------------
// Description
// This class calculate jet profiles for jets in a jet collection
// Author: L. Cunquiero, K. Garner, M. Verweij
//------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include "TF1.h"
#include "TMath.h"
#include "TH2.h"

#include "fastjet/PseudoJet.hh"
#include "jetCollection.hh"

class jetProfile {
private:

  std::vector<fastjet::PseudoJet> jets_;

  double eventWeight_;

  std::vector<double> min_delR_;
  std::vector<double> max_delR_;
  
  std::vector<std::vector<double>> jetProfiles_;
  
public:
  jetProfile(jetCollection &c, double eventWeight = 1.);
  jetProfile(std::vector<fastjet::PseudoJet> v, double eventWeight = 1.);

  TH2F* calculateProfileHisto();
  void  calculateProfile();

  void setBoundariesMin(std::vector<double> v) { min_delR_ = v; }
  void setBoundariesMax(std::vector<double> v) { max_delR_ = v; }
  
  std::vector<double> getJetProfile(int ijet) const {return jetProfiles_[ijet]; }
  std::vector<std::vector<double>> getJetProfiles() const {return jetProfiles_; }
  
};

  
jetProfile::jetProfile(jetCollection &c, double eventWeight) :
  jets_(c.getJet()),
  eventWeight_(eventWeight)
{
}

jetProfile::jetProfile(std::vector<fastjet::PseudoJet> v, double eventWeight) :
  jets_(v),
  eventWeight_(eventWeight)
{
}
  
TH2F* jetProfile::calculateProfileHisto()
{
  TH2F *h2PtJetRho = new TH2F("h2PtJetRhotemp","",150,0,150,20,0,0.5);
  double min_pPt = 1;
  double max_delR  = 0.5;

  for(fastjet::PseudoJet& jet : jets_) {
  
    std::vector<fastjet::PseudoJet> const_all = jet.constituents();
    std::vector<fastjet::PseudoJet> const_charged, const_neutral;
    SelectorIsCharged().sift(const_all, const_charged, const_neutral);
     
    for(fastjet::PseudoJet& p : const_charged) {
      double pPt  = p.perp();
      double delR   = p.delta_R(jet);
      if(pPt > min_pPt && delR < max_delR) {
	h2PtJetRho->Fill(jet.perp(),delR,eventWeight_*(pPt/jet.perp()));
      }
    }
  }
  return h2PtJetRho;
}

void jetProfile::calculateProfile()
{
  const int np = min_delR_.size();
 
  jetProfiles_.clear();
  jetProfiles_.reserve(jets_.size());

  for(fastjet::PseudoJet& jet : jets_) {

    std::vector<double> jetProfile(np);
    std::fill(jetProfile.begin(),jetProfile.end(),0.);

    std::vector<fastjet::PseudoJet> const_all = jet.constituents();
    for(fastjet::PseudoJet& p : const_all) {

      double prof = 0.;
      if(jet.perp()>0.) prof = p.perp()/jet.perp();

      double delR   = p.delta_R(jet);
      for(int i = 0; i<np; ++i) {
        if(delR>=min_delR_[i] && delR<max_delR_[i]) {
          jetProfile[i] += prof;
        }
      }
    }

    for(int i = 0; i<np; ++i) {
      double dr = max_delR_[i]-min_delR_[i];
      double prof = jetProfile[i];
      if(dr>0.) jetProfile[i] = prof/dr;
      else      jetProfile[i] = 0.;
    }
    jetProfiles_.push_back(jetProfile);
  }//jet loop
}

#endif
