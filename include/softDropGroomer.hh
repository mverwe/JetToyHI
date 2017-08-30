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

#include "jetCollection.hh"

//---------------------------------------------------------------
// Description
// This class runs SoftDrop grooming on a set of jets
// Author: M. Verweij, Y. Chen
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
  std::vector<double>             dr12_;       //distance between the two subjets
  std::vector<std::vector<fastjet::PseudoJet>> constituents_;
  std::vector<std::vector<fastjet::PseudoJet>> constituents1_;
  std::vector<std::vector<fastjet::PseudoJet>> constituents2_;

public :
  softDropGroomer(double zcut = 0.1, double beta = 0., double r0 = 0.4);
  void setZcut(double c);
  void setBeta(double b);
  void setR0(double r);
  void setInputJets(std::vector<fastjet::PseudoJet> v);
  std::vector<fastjet::PseudoJet> getGroomedJets() const;
  std::vector<double> getZgs() const;
  std::vector<int> getNDroppedSubjets() const;
  std::vector<double> getDR12() const;
  std::vector<fastjet::PseudoJet> doGrooming(jetCollection &c);
  std::vector<fastjet::PseudoJet> doGrooming(std::vector<fastjet::PseudoJet> v);
  std::vector<fastjet::PseudoJet> doGrooming();
  std::vector<std::vector<fastjet::PseudoJet>> getConstituents() {return constituents_;}
  std::vector<std::vector<fastjet::PseudoJet>> getConstituents1() {return constituents1_;}
  std::vector<std::vector<fastjet::PseudoJet>> getConstituents2() {return constituents2_;}
};

softDropGroomer::softDropGroomer(double zcut, double beta, double r0)
   : zcut_(zcut), beta_(beta), r0_(r0)
{
}

void softDropGroomer::setZcut(double c)
{
   zcut_ = c;
}

void softDropGroomer::setBeta(double b)
{
   beta_ = b;
}

void softDropGroomer::setR0(double r)
{
   r0_ = r;
}

void softDropGroomer::setInputJets(const std::vector<fastjet::PseudoJet> v)
{
   fjInputs_ = v;
}

std::vector<fastjet::PseudoJet> softDropGroomer::getGroomedJets() const
{
   return fjOutputs_;
}

std::vector<double> softDropGroomer::getZgs() const
{
   return zg_;
}

std::vector<int> softDropGroomer::getNDroppedSubjets() const
{
   return drBranches_;
}

std::vector<double> softDropGroomer::getDR12() const
{
   return dr12_;
}

std::vector<fastjet::PseudoJet> softDropGroomer::doGrooming(jetCollection &c)
{
   return doGrooming(c.getJet());
}

std::vector<fastjet::PseudoJet> softDropGroomer::doGrooming(std::vector<fastjet::PseudoJet> v)
{
   setInputJets(v);
   return doGrooming();
}

std::vector<fastjet::PseudoJet> softDropGroomer::doGrooming()
{
   fjOutputs_.reserve(fjInputs_.size());
   zg_.reserve(fjInputs_.size());
   drBranches_.reserve(fjInputs_.size());
   dr12_.reserve(fjInputs_.size());
   constituents_.reserve(fjInputs_.size());
   constituents1_.reserve(fjInputs_.size());
   constituents2_.reserve(fjInputs_.size());
   
   for(fastjet::PseudoJet& jet : fjInputs_) {
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

      fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 999.);
      fastjet::ClusterSequence cs(particles, jet_def);

      std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets());
      if(tempJets.size()<1) {
         fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         zg_.push_back(-1.);
         drBranches_.push_back(-1.);
         constituents1_.push_back(std::vector<fastjet::PseudoJet>());
         constituents2_.push_back(std::vector<fastjet::PseudoJet>());
         continue;
      }

      fastjet::contrib::SoftDrop * sd = new fastjet::contrib::SoftDrop(beta_, zcut_, r0_ );
      sd->set_verbose_structure(true);

      fastjet::PseudoJet transformedJet = tempJets[0];
      if ( transformedJet == 0 ) {
         fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         zg_.push_back(-1.);
         drBranches_.push_back(-1.);
         constituents1_.push_back(std::vector<fastjet::PseudoJet>());
         constituents2_.push_back(std::vector<fastjet::PseudoJet>());
         if(sd) { delete sd; sd = 0;}
         continue;
      }
      transformedJet = (*sd)(transformedJet);

      fastjet::PseudoJet j1, j2;
      if(transformedJet.has_parents(j1, j2) == true)
      {
         constituents1_.push_back(j1.constituents());
         constituents2_.push_back(j2.constituents());
      }
      else
      {
         constituents1_.push_back(std::vector<fastjet::PseudoJet>());
         constituents2_.push_back(std::vector<fastjet::PseudoJet>());
      }
      constituents_.push_back(transformedJet.constituents());

      double zg = transformedJet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
      int ndrop = transformedJet.structure_of<fastjet::contrib::SoftDrop>().dropped_count();
      //std::vector<double> dropped_symmetry = transformedJet.structure_of<fastjet::contrib::SoftDrop>().dropped_symmetry();
            
      //get distance between the two subjets
      std::vector<fastjet::PseudoJet> subjets;
      if ( transformedJet.has_pieces() ) {
        subjets = transformedJet.pieces();
        fastjet::PseudoJet subjet1 = subjets[0];
        fastjet::PseudoJet subjet2 = subjets[1];
        dr12_.push_back(subjet1.delta_R(subjet2));
      } else {
        dr12_.push_back(-1.);
      }

      fjOutputs_.push_back( transformedJet ); //put CA reclusterd jet after softDrop into vector
      zg_.push_back(zg);
      drBranches_.push_back(ndrop);

      if(sd) { delete sd; sd = 0;}
   }
   return fjOutputs_;
}

#endif
