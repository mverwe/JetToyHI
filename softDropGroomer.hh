#ifndef softDropGroomer_h
#define softDropGroomer_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "fastjet/contrib/SoftDrop.hh"

using namespace std;

//---------------------------------------------------------------
// Description
// This class runs SoftDrop grooming on a set of jets
// Author: M. Verweij
//---------------------------------------------------------------

class softDropGroomer {

private :
  double zcut_;
  double beta_;
  double r0_;
  
  std::vector<fastjet::PseudoJet> fjInputs_;   //ungroomed jets
  std::vector<fastjet::PseudoJet> fjOutputs_;  //groomed jets
  std::vector<double>             zg_;         //zg of groomed jets
  std::vector<int>                drBranches_; //dropped branches
  
public :
  softDropGroomer(double zcut = 0.1, double beta = 0., double r0 = 0.4) :
    zcut_(zcut), beta_(beta), r0_(r0)
  {
  }

  void setZcut(double c)     { zcut_ = c; }
  void setBeta(double b)     { beta_ = b; }
  void setR0(double r)       { r0_   = r; }

  void setInputJets(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }

  std::vector<fastjet::PseudoJet> getGroomedJets() const { return fjOutputs_;}
  std::vector<double> getZgs()                     const { return zg_;}
  std::vector<int> getNDroppedBranches()           const { return drBranches_;}
  
  std::vector<fastjet::PseudoJet> doGrooming(std::vector<fastjet::PseudoJet> v)
  {
     setInputJets(v);
     return doGrooming();
  }
  
  std::vector<fastjet::PseudoJet> doGrooming()
  {
    fjOutputs_.reserve(fjInputs_.size());
    zg_.reserve(fjInputs_.size());
    drBranches_.reserve(fjInputs_.size());
    
    for(fastjet::PseudoJet& jet : fjInputs_) {
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

      fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm,999.);
      fastjet::ClusterSequence cs(particles, jet_def);

      std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets());
      if(tempJets.size()<1) {
        fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
        zg_.push_back(-1.);
        drBranches_.push_back(-1.);
        continue;
      }
      
      fastjet::contrib::SoftDrop * sd = new fastjet::contrib::SoftDrop(beta_, zcut_, r0_ );
      sd->set_verbose_structure(true);

      fastjet::PseudoJet transformedJet = tempJets[0];
      if ( transformedJet == 0 ) {
        fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
        zg_.push_back(-1.);
        drBranches_.push_back(-1.);
        if(sd) { delete sd; sd = 0;}
        continue;
      }
      transformedJet = (*sd)(transformedJet);

      double zg = transformedJet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
      int ndrop = transformedJet.structure_of<fastjet::contrib::SoftDrop>().dropped_count();
      //std::vector<double> dropped_symmetry = transformedJet.structure_of<fastjet::contrib::SoftDrop>().dropped_symmetry();
      // std::cout << "ptG: " << transformedJet.pt() << " zg: " << zg << std::endl;
      
      fjOutputs_.push_back( transformedJet ); //put CA reclusterd jet after softDrop into vector
      zg_.push_back(zg);
      drBranches_.push_back(ndrop);

      if(sd) { delete sd; sd = 0;}
    }
    return fjOutputs_;
  }
};

#endif
