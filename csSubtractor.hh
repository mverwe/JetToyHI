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

using namespace std;

//---------------------------------------------------------------
// Description
// This class runs the jet-by-jet constituent subtraction
// Author: M. Verweij
//---------------------------------------------------------------

class csSubtractor {

private :
  double alpha_;
  double rParam_;
  double ghostArea_;
  double ghostRapMax_;
  double jetRapMax_;
  double rho_;
  double rhom_;
  std::vector<fastjet::PseudoJet> fjInputs_;

  contrib::ConstituentSubtractor subtractor_;

  
public :
  csSubtractor(double alpha = 1., double rParam = -1., double ghostArea = 0.005, double ghostRapMax = 3.0, double jetRapMax = 3.0) :
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
    subtractor_.set_do_mass_subtraction(true);

  }

  void setAlpha(double a) { alpha_ = a; }
  void setRParam(double r) { rParam_ = r; }
  void setGhostArea(double a) { ghostArea_ = a; }

  void setInputParticles(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }

  double getRho()  const { return rho_; }
  double getRhoM() const { return rhom_; }
  
  std::vector<fastjet::PseudoJet> doSubtraction() {

    // do the clustering with ghosts and get the jets
    //----------------------------------------------------------
    fastjet::JetDefinition jet_def(antikt_algorithm, rParam_);
    fastjet::GhostedAreaSpec ghost_spec(ghostRapMax_, 1, ghostArea_);
    fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,ghost_spec);
    
    fastjet::ClusterSequenceArea cs(fjInputs_, jet_def, area_def);
    fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax_);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(jet_selector(cs.inclusive_jets()));
    
    // create what we need for the background estimation
    //----------------------------------------------------------
    fastjet::JetDefinition jet_def_bkgd(fastjet::kt_algorithm, 0.4);
    fastjet::AreaDefinition area_def_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
    fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4) * (!fastjet::SelectorNHardest(2));
    fastjet::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
    bkgd_estimator.set_particles(fjInputs_);

    rho_ = bkgd_estimator.rho();
    rhom_ = bkgd_estimator.rho_m();

    /*
    cout << "Background estimation:" << endl;
    cout << "  " << bkgd_estimator.description() << endl << endl;;
    cout << "  Giving, for the full event" << endl;
    cout << "    rho     = " << bkgd_estimator.rho()   << endl;
    cout << "    sigma   = " << bkgd_estimator.sigma() << endl;
#if FASTJET_VERSION_NUMBER >= 30100
    cout << "    rho_m   = " << bkgd_estimator.rho_m()   << endl;
    cout << "    sigma_m = " << bkgd_estimator.sigma_m() << endl;
#endif
    cout << endl;
    */
    subtractor_.set_background_estimator(&bkgd_estimator);
    subtractor_.set_common_bge_for_rho_and_rhom(true);

    std::vector<fastjet::PseudoJet> csjets;
    csjets.reserve(jets.size());
    for(fastjet::PseudoJet& jet : jets) {
      fastjet::PseudoJet subtracted_jet = subtractor_(jet);
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(subtracted_jet.constituents(), ghosts, particles);
      if(particles.size()>0.) //don't store ghost jets
        csjets.push_back(subtracted_jet);
    }

    return csjets;
  }
};

#endif
