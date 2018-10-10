#ifndef sharedLayerSubtractor_h
#define sharedLayerSubtractor_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <map>
#include <random>
#include <numeric>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "../PU14/PU14.hh"

#include "Angularity.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs the jet-by-jet shared layer subtraction
// Author: M. Verweij, Y. Mehtar-Tani
//---------------------------------------------------------------

class sharedLayerSubtractor {

private :
  double jetRParam_;
  double ghostArea_;
  double ghostRapMax_;
  double jetRapMax_;
  double rho_;
  double rhoSigma_;
  double pTDbkg_;
  double pTDbkgSigma_;
  
  std::vector<fastjet::PseudoJet> fjInputs_;
  std::vector<fastjet::PseudoJet> fjJetInputs_;

  std::vector<fastjet::PseudoJet> fjJetParticles_;

  std::random_device rd_;
  int nInitCond_;
  int nTopInit_;
    
public :
  sharedLayerSubtractor(double rJet = 0.4,
                        double ghostArea = 0.005,
                        double ghostRapMax = 3.0,
                        double jetRapMax = 3.0,
                        int nInitCond = 100,
                        int nTopInit = 10) :
    jetRParam_(rJet),
    ghostArea_(ghostArea),
    ghostRapMax_(ghostRapMax),
    jetRapMax_(jetRapMax),
    nInitCond_(nInitCond),
    nTopInit_(nTopInit)
  {

  }

  void setGhostArea(double a) { ghostArea_ = a; }

  void setInputParticles(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }
  void setInputJets(std::vector<fastjet::PseudoJet> v)      { fjJetInputs_ = v; }

  double getRho()  const { return rho_; }
  double getRhoSigma() const { return rhoSigma_; }

  double getPTDBkg() const { return pTDbkg_; }
  double getPTDBkgSigma() const { return pTDbkgSigma_; }

  std::vector<fastjet::PseudoJet> doSubtraction() {

    //if(fjJetInputs_.size()==0 && fjInputs_.size()) {
    //  throw "You didn't give me input jets or particles. You should give me one of the two";
    //  return std::vector<fastjet::PseudoJet>();
    //}

    fastjet::GhostedAreaSpec ghost_spec(ghostRapMax_, 1, ghostArea_);

    std::vector<fastjet::PseudoJet> jets = fjJetInputs_;
    //if(jets.size()==0) {

    // do the clustering with ghosts and get the jets
    //----------------------------------------------------------
    fastjet::JetDefinition jet_def(antikt_algorithm, jetRParam_);
    fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,ghost_spec);

    fastjet::ClusterSequenceArea cs(fjInputs_, jet_def, area_def);
    fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax_);
    jets = fastjet::sorted_by_pt(jet_selector(cs.inclusive_jets()));

    fastjet::JetDefinition jet_defSub(antikt_algorithm, 999.);
    fastjet::ClusterSequenceArea csSub(fjInputs_, jet_defSub, area_def);
    
    // create what we need for the background estimation
    //----------------------------------------------------------
    fastjet::JetDefinition jet_def_bkgd(fastjet::kt_algorithm, 0.4);
    fastjet::AreaDefinition area_def_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
    fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4) * (!fastjet::SelectorNHardest(2));

    fastjet::ClusterSequenceArea csKt(fjInputs_, jet_def_bkgd, area_def_bkgd);
    std::vector<fastjet::PseudoJet> bkgd_jets = fastjet::sorted_by_pt(selector(csKt.inclusive_jets()));
    
    fastjet::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
    bkgd_estimator.set_particles(fjInputs_);
    bkgd_estimator.set_jets(bkgd_jets);
    
    rho_ = bkgd_estimator.rho();
    rhoSigma_ = bkgd_estimator.sigma();

    if(rho_ < 0)    rho_ = 0;

    //std::cout << "rho: " << rho_ << "  rhoSigma: " << rhoSigma_ << std::endl;
    
    //initial gaus with mean=rho_ and width=rhoSigma_
    std::normal_distribution<> gausDist(rho_, rhoSigma_);

    //UE metric
    //----------------------------------------------------------
    //Angularity width(1.,1.,0.4);
    Angularity pTD(0.,2.,0.4);

    std::vector<double> pTD_bkgd;
    for(fastjet::PseudoJet& jet : bkgd_jets) {
      pTD_bkgd.push_back(pTD.result(jet));
    }
    std::nth_element(pTD_bkgd.begin(), pTD_bkgd.begin() + pTD_bkgd.size()/2, pTD_bkgd.end());
    double med_pTD = pTD_bkgd[pTD_bkgd.size()/2];

    int nRMS = 0;
    double rms_pTD = 0.;
    for(int ip = 0; ip<(int)pTD_bkgd.size(); ++ip) {
      rms_pTD += (pTD_bkgd[ip]-med_pTD)*(pTD_bkgd[ip]-med_pTD);
      nRMS++;
    }
    if(nRMS>0.)
      rms_pTD = sqrt(rms_pTD/(double)nRMS);

    pTDbkg_ = med_pTD;
    pTDbkgSigma_ = rms_pTD;

    //std::cout << "med_pTD: " << med_pTD << "  rms_pTD: " << rms_pTD << std::endl;
    
    std::mt19937 rndSeed(rd_()); //rnd number generator seed
           
    std::vector<fastjet::PseudoJet> subtracted_jets;
    subtracted_jets.reserve(jets.size());
    int ijet = -1;
    for(fastjet::PseudoJet& jet : jets) {
      ++ijet;
      //if(ijet>1) continue; //TEMP,DEBUG
      // std::cout << "start jet loop. entry: " << ijet << "/" << jets.size() << " pt: " << jet.pt() << " eta: " << jet.eta() << std::endl;

      if(jet.is_pure_ghost()) continue;
      
      //get ghosts and true particles (ghosts are distributed uniformly which we will use to create initial conditions)
      //----------------------------------------------------------
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
      if(particles.size()<1 || jet.pt()<1.) continue;
      
      //set user_index of all particles to position particles vector
      for(int i = 0; i<(int)particles.size(); ++i) {
        particles[i].set_user_index(i);
      }
      
      // create requested number of initial conditions
      //----------------------------------------------------------
      std::vector<std::vector<int>> collInitCond;
      for(int ii = 0; ii<nInitCond_; ++ii) {

        std::uniform_int_distribution<> distUni(0,ghosts.size()); //uniform distribution of ghosts in vector
        std::vector<int> initCondition;                           //list of particles in initial condition

        //get random maxPt for this initial condition
        double maxPt = gausDist(rndSeed)*jet.area();

        //make copy of particles so that a particle is not repeated inside the same initial condition
        std::vector<fastjet::PseudoJet> particlesNotUsed = particles;
        
        double maxPtCurrent = 0.;
        while(maxPtCurrent<maxPt && particlesNotUsed.size()>0) {

          //pick random ghost
          int ighost = int(std::floor(distUni(rndSeed))); 
        
          //find closest particle to ghost
          double dR2 = 1000.;
          fastjet::PseudoJet partSel;
          int ipSel = -1;
          for(int ip = 0; ip<(int)particlesNotUsed.size(); ++ip) {
            fastjet::PseudoJet part = particlesNotUsed[ip];
            double dR2Cur = ghosts[ighost].squared_distance(part);
            if(dR2Cur < dR2) {
              dR2 = dR2Cur;
              partSel = part;
              ipSel = ip;
            }
          }
          if(ipSel<0) continue; //this shouldn't happen
          initCondition.push_back(partSel.user_index());
          maxPtCurrent+=partSel.pt();
          particlesNotUsed.erase(particlesNotUsed.begin()+ipSel);
          //std::cout << "Added new particle with pt = " << partSel.pt() << " to init condition. total pt now " << maxPtCurrent << "/" << maxPt << std::endl;
        }
        collInitCond.push_back(initCondition);
      }//initial conditions loop
     
      //----------------------------------------------------------
      //Now we have the requested number of random initial condition
      
      //Next step: calc chi2 for each initial condition
      //----------------------------------------------------------
      std::vector<double> chi2s;
      for(int ii = 0; ii<nInitCond_; ++ii) {
        std::vector<int> indices = collInitCond[ii];
        double chi2 = 1e6;
        if(indices.size()>0) { 
          std::vector<fastjet::PseudoJet> combinedparticles;
          for(int ic = 0; ic<(int)indices.size(); ++ic) {
            combinedparticles.push_back(particles[indices[ic]]);
          }
          fastjet::PseudoJet currInitJet = fastjet::PseudoJet(join(combinedparticles));
          //std::cout << "get pTD from ID" << std::endl;
          double ptDCur = pTD.result(currInitJet);
          chi2 = fabs(ptDCur-med_pTD)*(fabs(ptDCur-med_pTD))/rms_pTD/rms_pTD;
        }
        chi2s.push_back(chi2);
        //std::cout << "chi2 of initial condition " << ii << " is " << chi2 << " with " <<  combinedparticles.size() << " particles" << std::endl;
      }
      
      
      //sort the chi2s keeping track of indices
      //----------------------------------------------------------
      // initialize original index locations
      std::vector<size_t> idx(chi2s.size());
      iota(idx.begin(), idx.end(), 0);

      // sort indexes based on comparing values in v
      std::sort(idx.begin(), idx.end(),
                [&chi2s](size_t i1, size_t i2) {return chi2s[i1] < chi2s[i2];});

            
      //Next step: figure out how often each particle is shared in nTopInit_ initial conditions
      //----------------------------------------------------------
      std::vector<int> share_idx(particles.size(),0);
      for(int it = 0; it<nTopInit_; ++it) {
        int chi2Index = idx[it];
        std::vector<int> indices = collInitCond[chi2Index];
        for(int ic = 0; ic<(int)indices.size(); ++ic) {
          share_idx[particles[indices[ic]].user_index()]++;
          //std::cout << "indices[ic] = " << indices[ic] << "  user_index = " << particles[indices[ic]].user_index() << std::endl;
          //std::cout << "share_idx current: " << share_idx[particles[indices[ic]].user_index()] << std::endl;
        }
      }
      
      //sort according to how often a particle is shared
      //----------------------------------------------------------
      // initialize original index locations
      std::vector<size_t> ish(share_idx.size());
      iota(ish.begin(), ish.end(), 0);
      std::sort(ish.begin(), ish.end(),
                [&share_idx](size_t i1, size_t i2) {return share_idx[i1] > share_idx[i2];});

      
      //create final UE object
      //----------------------------------------------------------
      double maxPtFinalUE = gausDist(rndSeed)*jet.area();
      double curPtFinalUE = 0.;
      std::vector<fastjet::PseudoJet> bkgd_particles;
      fjJetParticles_.clear();            
      for(auto userIndex : ish) {
        //std::cout << "userIndex: " << userIndex << " nshared: " << share_idx[userIndex] << std::endl;
        fastjet::PseudoJet part = particles[userIndex]; 
        if(curPtFinalUE<maxPtFinalUE) { //assign as bkgd particle
          curPtFinalUE+=part.pt();
          bkgd_particles.push_back(part);
        } else { //assign as jet particle
          fjJetParticles_.push_back(part);
        }
      }

      if(fjJetParticles_.size()>0) {
        fastjet::ClusterSequenceArea *csSub = new fastjet::ClusterSequenceArea(fjJetParticles_, jet_defSub, area_def);
        std::vector<fastjet::PseudoJet> jetSub = fastjet::sorted_by_pt(csSub->inclusive_jets());
        if(jetSub[0].pt()>0.) subtracted_jets.push_back(jetSub[0]);
        if(subtracted_jets.size()>0 && subtracted_jets.size()<2) csSub->delete_self_when_unused();
      }
    }//jet loop

    //std::cout << "\n n subtracted jets: " << subtracted_jets.size() << std::endl;
    return subtracted_jets;
  }
};

#endif
