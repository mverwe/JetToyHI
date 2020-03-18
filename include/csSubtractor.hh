#ifndef csSubtractor_h
#define csSubtractor_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "fastjet/contrib/ConstituentSubtractor.hh"

#include "../PU14/PU14.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs the jet-by-jet constituent subtraction
// Author: M. Verweij, Y. Chen
//---------------------------------------------------------------

class csSubtractor {

   private :
      double jetRParam_;
      double alpha_;
      double rParam_;
      double ghostArea_;
      double ghostRapMax_;
      double jetRapMax_;
      double rho_;
      double rhom_;
      std::vector<fastjet::PseudoJet> fjInputs_;
      std::vector<fastjet::PseudoJet> fjJetInputs_;
      std::vector<std::vector<fastjet::PseudoJet>> Hard;
      std::vector<std::vector<fastjet::PseudoJet>> Soft;

      contrib::ConstituentSubtractor subtractor_;


   public :
      csSubtractor(double rJet = 0.4, double alpha = 1., double rParam = -1., double ghostArea = 0.005, double ghostRapMax = 3.0, double jetRapMax = 3.0) :
         jetRParam_(rJet),
         alpha_(alpha),
         rParam_(rParam),
         ghostArea_(ghostArea),
         ghostRapMax_(ghostRapMax),
         jetRapMax_(jetRapMax)
   {
      //init constituent subtractor
      subtractor_.set_distance_type(contrib::ConstituentSubtractor::deltaR);
      subtractor_.set_max_distance(rParam_); //free parameter for the maximal allowed distance between particle i and ghost k
      subtractor_.set_alpha(alpha_); // free parameter for the distance measure (the exponent of particle pt). Note that in older versions of the package alpha was multiplied by two but in newer versions this is not the case anymore
      subtractor_.set_do_mass_subtraction();//true);

   }

      void setAlpha(double a) { alpha_ = a; }
      void setRParam(double r) { rParam_ = r; }
      void setGhostArea(double a) { ghostArea_ = a; }

      void setInputParticles(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }
      void setInputJets(std::vector<fastjet::PseudoJet> v)      { fjJetInputs_ = v; }

      double getRho()  const { return rho_; }
      double getRhoM() const { return rhom_; }
      std::vector<std::vector<fastjet::PseudoJet>> getHard() const { return Hard; }
      std::vector<std::vector<fastjet::PseudoJet>> getSoft() const { return Soft; }

      std::vector<fastjet::PseudoJet>  getUnsubtractedJets() const { return fjJetInputs_;}

      std::vector<fastjet::PseudoJet> doSubtraction() {

         Hard.clear();
         Soft.clear();

         //if(fjJetInputs_.size()==0 && fjInputs_.size()) {
         //  throw "You didn't give me input jets or particles. You should give me one of the two";
         //  return std::vector<fastjet::PseudoJet>();
         //}

 
         std::vector<fastjet::PseudoJet> jets = fjJetInputs_;
         //if(jets.size()==0) {

         // do the clustering with ghosts and get the jets
         //----------------------------------------------------------
         fastjet::GhostedAreaSpec ghost_spec(ghostRapMax_, 1, ghostArea_);
         fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,ghost_spec);
         fastjet::JetDefinition jet_def(antikt_algorithm, jetRParam_);
                           
         fastjet::ClusterSequenceArea cs(fjInputs_, jet_def, area_def);
         fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax_);
         jets = fastjet::sorted_by_pt(jet_selector(cs.inclusive_jets()));
         // fjJetInputs_ = jets;
         //}

         // create what we need for the background estimation
         //----------------------------------------------------------
         fastjet::JetDefinition jet_def_bkgd(fastjet::kt_algorithm, 0.4);
         fastjet::AreaDefinition area_def_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
         fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4) * (!fastjet::SelectorNHardest(2));
         fastjet::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
         bkgd_estimator.set_particles(fjInputs_);

         rho_ = bkgd_estimator.rho();
         rhom_ = bkgd_estimator.rho_m();

         if(rho_ < 0)    rho_ = 0;
         if(rhom_ < 0)   rhom_ = 0;

         subtractor_.set_background_estimator(&bkgd_estimator);
         subtractor_.set_common_bge_for_rho_and_rhom(true);

         std::vector<fastjet::PseudoJet> csjets;
         csjets.reserve(jets.size());
         for(fastjet::PseudoJet& jet : jets) {

            fastjet::PseudoJet subtracted_jet = subtractor_(jet);
            std::vector<fastjet::PseudoJet> particles, ghosts;
            fastjet::SelectorIsPureGhost().sift(subtracted_jet.constituents(), ghosts, particles);

            std::vector<fastjet::PseudoJet> A, B;
            SelectorIsHard().sift(jet.constituents(), A, B);

            std::vector<fastjet::PseudoJet> hard, soft;
            for(fastjet::PseudoJet p : subtracted_jet.constituents())
            {
               double BestA = -1, BestB = -1;
               for(fastjet::PseudoJet x : A)
               {
                  double DR2 = p.squared_distance(x);
                  if(BestA < 0 || BestA > DR2)
                     BestA = DR2;
               }
               for(fastjet::PseudoJet x : B)
               {
                  double DR2 = p.squared_distance(x);
                  if(BestB < 0 || BestB > DR2)
                     BestB = DR2;
               }
               if(BestA < BestB)
                  hard.push_back(p);
               else
                  soft.push_back(p);
            }

            if(particles.size() > 0)
            {
               std::vector<fastjet::PseudoJet> combinedparticles;
               for(fastjet::PseudoJet p : particles)
                  combinedparticles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.E()));
               for(fastjet::PseudoJet p : jet.constituents())
                  if(p.E() < 1e-5)
                     combinedparticles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.E()));
         
               csjets.push_back(fastjet::PseudoJet(join(combinedparticles)));
               Hard.push_back(hard);
               Soft.push_back(soft);
            }
         }

         return csjets;
      }
};

#endif
