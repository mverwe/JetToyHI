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

class softDropGroomer {

private :
  double zcut_;
  double beta_;
  double r0_;
  double jetPtMin_;
  
  std::vector<fastjet::PseudoJet> fjInputs_;   //ungroomed jets
  std::vector<fastjet::PseudoJet> fjOutputs_;  //groomed jets
  std::vector<double>             zg_;         //zg of groomed jets
  std::vector<int>                drBranches_; //dropped branches
  
public :
  softDropGroomer() {
    zcut_ = 0.1;
    beta_ = 0.;
    r0_   = 0.4;
    jetPtMin_ = 10.;
  }

  void setZcut(double c)     { zcut_ = c; }
  void setBeta(double b)     { beta_ = b; }
  void setR0(double r)       { r0_   = r; }
  void setJetPtMin(double m) { jetPtMin_ = m; }

  void setInputJets(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }

  std::vector<fastjet::PseudoJet> getGroomedJets() const { return fjOutputs_;}
  std::vector<double> getZgs()                     const { return zg_;}
  std::vector<int> getNDroppedBranches()           const { return drBranches_;}
  
  std::vector<fastjet::PseudoJet> doGrooming() {

    fjOutputs_.reserve(fjInputs_.size());
    zg_.reserve(fjInputs_.size());
    drBranches_.reserve(fjInputs_.size());
    for(fastjet::PseudoJet& jet : fjInputs_) {
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

      fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm,999.);
      fastjet::ClusterSequence cs(particles, jet_def);

      std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets(jetPtMin_));
      if(tempJets.size()<1) continue;

      fastjet::contrib::SoftDrop * sd = new fastjet::contrib::SoftDrop(beta_, zcut_, r0_ );
      sd->set_verbose_structure(true);

      fastjet::PseudoJet transformedJet = tempJets[0];
      if ( transformedJet == 0 ) {
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
