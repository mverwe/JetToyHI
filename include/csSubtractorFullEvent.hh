#ifndef csSubtractorFullEvent_h
#define csSubtractorFullEvent_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "fastjet/contrib/ConstituentSubtractor.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs the full event constituent subtraction
// Author: M. Verweij
//---------------------------------------------------------------

class csSubtractorFullEvent {

private :
  double alpha_;
  double rParam_;
  double ghostArea_;
  double ghostRapMax_;
  double rho_;
  double rhom_;
  std::vector<fastjet::PseudoJet> fjInputs_;

  contrib::ConstituentSubtractor subtractor_;

  
public :
  csSubtractorFullEvent(double alpha = 1., double rParam = 0.4, double ghostArea = 0.005, double ghostRapMax = 3.0) :
    alpha_(alpha),
    rParam_(rParam),
    ghostArea_(ghostArea),
    ghostRapMax_(ghostRapMax),
    rho_(-1),
    rhom_(-1)
  {
    //init constituent subtractor
    subtractor_.set_distance_type(contrib::ConstituentSubtractor::deltaR);
    subtractor_.set_max_distance(rParam_); //free parameter for the maximal allowed distance between particle i and ghost k
    subtractor_.set_alpha(alpha_); // free parameter for the distance measure (the exponent of particle pt). Note that in older versions of the package alpha was multiplied by two but in newer versions this is not the case anymore
    subtractor_.set_scale_fourmomentum(); //Keep rapidity and pseudo-rapidity fixed (scale fourmomentum). Recommended - observed better performance than the mass correction. Use: subtractor.set_scale_fourmomentum();
    subtractor_.set_do_mass_subtraction();//true);

  }

  void setAlpha(double a)     { alpha_ = a; }
  void setRParam(double r)    { rParam_ = r; }
  void setGhostArea(double a) { ghostArea_ = a; }

  void setRho(double r)       { rho_ = r; }
  void setRhom(double r)      { rhom_ = r; }

  void setInputParticles(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }

  double getRho()  const { return rho_; }
  double getRhoM() const { return rhom_; }
  
  std::vector<fastjet::PseudoJet> doSubtraction() {

    if(rho_<0.) {    
      // create what we need for the background estimation
      //----------------------------------------------------------
      fastjet::GhostedAreaSpec ghost_spec(ghostRapMax_, 1, ghostArea_);
      fastjet::JetDefinition jet_def_bkgd(fastjet::kt_algorithm, 0.4);
      fastjet::AreaDefinition area_def_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
      fastjet::Selector selector = fastjet::SelectorAbsRapMax(ghostRapMax_-0.4) * (!fastjet::SelectorNHardest(2));
      fastjet::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
      bkgd_estimator.set_particles(fjInputs_);
      
      rho_ = bkgd_estimator.rho();
      rhom_ = bkgd_estimator.rho_m();
      
      subtractor_.set_background_estimator(&bkgd_estimator);
      subtractor_.set_common_bge_for_rho_and_rhom(true);
    } else {
      //if rho and rhom provided, use externally supplied densities
      subtractor_ = contrib::ConstituentSubtractor(rho_,rhom_,alpha_,rParam_,contrib::ConstituentSubtractor::deltaR);
    }
    
    std::vector<fastjet::PseudoJet> corrected_event = subtractor_.subtract_event(fjInputs_,ghostRapMax_);
    
    return corrected_event;
  }
};

#endif
