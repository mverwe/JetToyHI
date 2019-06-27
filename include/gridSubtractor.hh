#ifndef gridSubtractor_h
#define gridSubtractor_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "jetCollection.hh"
#include "jewelMatcher.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs the GridSub1 method for JEWEL with recoils following
// httphs:jewel.hepforge.org/JEWEL_BKGSUBTRACTION.cc
// Author: A. Soto-Ontoso
//---------------------------------------------------------------

class gridSubtractor {

private :
  std::vector<fastjet::PseudoJet> fjScatteringCenters_; // Scattering Centers
  std::vector<fastjet::PseudoJet> fjInputs_;   //unsubtracted jets
  std::vector<fastjet::PseudoJet> fjOutputs_;  //subtracted jets

  std::vector<fastjet::PseudoJet> pJet_sub; // list of particles after subtraction (this is the output)

  double etaMin_;
  double phiMin_;
  double etaMax_;
  double phiMax_;
  double jetRParam_;
  double deltaRMin_;

public :
  gridSubtractor(double etaMin = -1., double phiMin =-M_PI, double etaMax = 1, double phiMax = M_PI, double jetRParam = 0.4, double deltaRMin = 0.05);
  void setInputJets(std::vector<fastjet::PseudoJet> v);
  void setParticlesCenters(std::vector<fastjet::PseudoJet> ScatteringCenters);
  std::vector<fastjet::PseudoJet> doGridSub1(jetCollection c, std::vector<fastjet::PseudoJet> ScatteringCenters);
  std::vector<fastjet::PseudoJet> doGridSub1();
  std::pair<int,int> FindCandidateBin(double eta, double phi, std::vector<double> _etabins, vector<double> _phibins);

  struct candidate{
    //! these are the boundaries
    double etaMin;
    double phiMin;
    double etaMax;
    double phiMax;
    //! Pseudo jet will have the objects inside
    std::vector<fastjet::PseudoJet> objects;
    std::vector<int> jetID;
    //! Four momentum to get useful informations like eta, phi, mass, pT etc...
    std::valarray<double> candMom = {0,0,0,0};
    std::valarray<double> bkgMom = {0,0,0,0};
    double sumnegpT = 0;
    double sumpT = 0;
  };
};

gridSubtractor::gridSubtractor(double etaMin, double phiMin, double etaMax, double phiMax, double jetRParam, double deltaRMin) : etaMin_(etaMin), phiMin_(phiMin), etaMax_(etaMax),phiMax_(phiMax),jetRParam_(jetRParam), deltaRMin_(deltaRMin)
{
}


void gridSubtractor::setInputJets(const std::vector<fastjet::PseudoJet> v)
{
  fjInputs_ = v;
}

void gridSubtractor::setParticlesCenters(const std::vector<fastjet::PseudoJet> ScatteringCenters)
{
  fjScatteringCenters_ = ScatteringCenters;
}

std::pair<int,int> gridSubtractor::FindCandidateBin(double eta, double phi, vector<double> _etabins, vector<double> _phibins){
  int etaloc = -9;
  int philoc = -9;

  for(unsigned x = 0; x<_etabins.size()-1; ++x){
    for(unsigned y = 0; y<_phibins.size()-1; ++y){
      if(eta >= _etabins[x] && eta < _etabins[x+1] &&
         phi >= _phibins[y] && phi < _phibins[y+1]){
        etaloc = x;
        philoc = y;
      }
    }
  }
  if(etaloc == -9 || philoc == -9)
    return std::pair<int,int>(-1,-1);
  else
    return std::pair<int,int>(etaloc,philoc);
}


std::vector<fastjet::PseudoJet> gridSubtractor::doGridSub1(jetCollection c, std::vector<fastjet::PseudoJet> ScatteringCenters)
{
  setInputJets(c.getJet());
  setParticlesCenters(ScatteringCenters);
  return doGridSub1();
}

std::vector<fastjet::PseudoJet> gridSubtractor::doGridSub1() {
  fjOutputs_.reserve(fjInputs_.size());
        
  //----------------------------------------------------
  // Create the grid (this could be somewhere else)
  //----------------------------------------------------

  // Cell boundaries
  std::vector<double> _etabins;
  std::vector<double> _phibins;

  int _Nbounds_eta = (int)2*etaMax_/(deltaRMin_);
  int _Nbounds_phi = (int)2*phiMax_/(deltaRMin_);

  for(int i = 0; i<=_Nbounds_eta; i++){
    double etax = -1*etaMax_ + i*deltaRMin_;
    _etabins.push_back(etax);
  }

  for(int i = 0; i<=_Nbounds_phi; i++){
    double phix = -1*phiMax_ + i*deltaRMin_;
    _phibins.push_back(phix);
  }

  //  std::cout<<"Grid we are using is ("<<_Nbounds_eta<<" x "<<_Nbounds_phi
  //		 <<")   eta: [-"<<etaMax_<<", "<<etaMax_
  //		 <<"]   phi: [-"<<phiMax_<<", "<<phiMax_<<"]"<<std::endl;

  // Initialize the grid

  std::vector<std::vector<candidate>> grid;
  for(int x = 0; x<_Nbounds_eta; x++){
    std::vector<candidate> xgrid;
    for(int y = 0; y<_Nbounds_phi; y++){
      candidate test;
      test.etaMin = _etabins[x];
      test.phiMin = _phibins[y];
      test.etaMax = _etabins[x+1];
      test.phiMax = _phibins[y+1];
      test.sumnegpT = 0.0;
      xgrid.push_back(test);
    }
    grid.push_back(xgrid);
    xgrid.clear();
  }

  //-----------------------------------------------------------------
  // Loop over all the jets and fill their constituents into the grid
  //-----------------------------------------------------------------

  int ijet = 0;

  for(fastjet::PseudoJet& jet : fjInputs_) {
    ijet++;
    std::vector<fastjet::PseudoJet> particles, ghosts;
    fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
    // loop over the candidates of tje jet and include them in the grid
    for (fastjet::PseudoJet p : particles){

      double peta = p.eta();
      double pphi = p.phi_std();
      std::valarray<double> p_fourmomentum = {p.px(), p.py(), p.pz(), p.E()};
      std::pair<int,int> candpos = FindCandidateBin(peta,pphi,_etabins,_phibins);
      if(candpos.first!=1 && candpos.second!=-1){
        grid[candpos.first][candpos.second].objects.push_back(p);
        grid[candpos.first][candpos.second].candMom+=p_fourmomentum;
        grid[candpos.first][candpos.second].jetID.push_back(ijet);
        grid[candpos.first][candpos.second].sumpT+=p.pt();
      }
    } // loop over constituents
  } // loop over jets

  //----------------------------------------------------------------
  // Loop over the scattering centers and add them as background
  //----------------------------------------------------------------

  for(fastjet::PseudoJet p : fjScatteringCenters_) {
    bool iJ = false;

    for(fastjet::PseudoJet& jet : fjInputs_) {
      double delR = p.delta_R(jet);
      if(delR < jetRParam_){
        iJ = true;
        break;
      }
    } //loop over jets

    if (iJ){
      std::valarray<double> p_fourmomentum = {p.px(), p.py(), p.pz(), p.E()};
      std::pair<int,int> candpos = FindCandidateBin(p.eta(),p.phi_std(),_etabins,_phibins);
      if(candpos.first!=1 && candpos.second!=-1){
        grid[candpos.first][candpos.second].objects.push_back(p);
        grid[candpos.first][candpos.second].bkgMom+=p_fourmomentum;
        grid[candpos.first][candpos.second].sumpT-=p.pt();
      }
    }

  } // loop over scattering centers

  //----------------------------------------------------------------
  // Build up the pseudojets to do the clustering
  //----------------------------------------------------------------

  std::vector<fastjet::PseudoJet> pJet_sub;

  for(int x = 0; x<_Nbounds_eta; x++){
    for (int y = 0; y<_Nbounds_phi; y++){
      std::valarray<double> jet_sub = grid[x][y].candMom;
      jet_sub-=grid[x][y].bkgMom;

      double px,py,pz,E = 0;

      double cand_px = grid[x][y].candMom[0];
      double cand_py = grid[x][y].candMom[1];
      double cand_pT = sqrt(pow(cand_px,2)+pow(cand_py,2));

      double bkg_px = grid[x][y].bkgMom[0];
      double bkg_py = grid[x][y].bkgMom[1];
      double bkg_pT = sqrt(pow(bkg_px,2)+pow(bkg_py,2));

      if((cand_pT-bkg_pT) < 0.0){
        grid[x][y].sumnegpT+=bkg_pT - cand_pT;
        px = 0;
        py = 0;
        pz = 0;
        E = 0;
      }
      else{
        px = jet_sub[0];
        py = jet_sub[1];
        pz = jet_sub[2];
        E = jet_sub[3];
      }
      PseudoJet part_sub(px,py,pz,E);
      if (part_sub.pt()!=0 && part_sub.E()!=0) pJet_sub.push_back(part_sub);
    } // y loop
  } // x loop


  return pJet_sub;
}


#endif
