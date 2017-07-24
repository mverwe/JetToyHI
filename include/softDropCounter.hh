#ifndef softDropCounter_h
#define softDropCounter_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Recluster.hh"

//---------------------------------------------------------------
// Description
// This class runs iterative SoftDrop on a set of jets
// Author: Y. Chen
//---------------------------------------------------------------

class softDropCounter
{
private :
   double zcut_;   // ZCut used in SD
   double beta_;   // Beta parameter
   double r0_;     // Jet radius
   double rcut_;   // Termination Criteria

   std::vector<fastjet::PseudoJet> fjInputs_;     // ungroomed jets
   std::vector<std::vector<double>> zgs_;         // all the zg's in the algorithm

public :
   softDropCounter(double z = 0.1, double beta = 0.0, double r0 = 0.4, double rcut = 0.1);
   void setZCut(double c);
   void setBeta(double b);
   void setR0(double r) ;
   void setRCut(double r);
   void setInputJets(const std::vector<fastjet::PseudoJet> &v);
   void run(const jetCollection &c);
   void run(const std::vector<fastjet::PseudoJet> &v);
   void run();
   std::vector<double> calculateNSD(double Kappa);
};

softDropCounter::softDropCounter(double z, double beta, double r0, double rcut)
   : zcut_(z), beta_(beta), r0_(r0), rcut_(rcut)
{
}

void softDropCounter::setZCut(double c)
{
   zcut_ = c;
}

void softDropCounter::setBeta(double b)
{
   beta_ = b;
}

void softDropCounter::setR0(double r)
{
   r0_   = r;
}

void softDropCounter::setRCut(double r)
{
   rcut_ = r;
}

void softDropCounter::setInputJets(const std::vector<fastjet::PseudoJet> &v)
{
   fjInputs_ = v;
}

void softDropCounter::run(const jetCollection &c)
{
   run(c.getJet());
}

void softDropCounter::run(const std::vector<fastjet::PseudoJet> &v)
{
   setInputJets(v);
   run();
}

void softDropCounter::run()
{
   int N = fjInputs_.size();

   for(fastjet::PseudoJet &jet: fjInputs_)
   {
      if(jet.has_constituents() == false)
      {
         zgs_.push_back(vector<double>());
         continue;
      }

      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

      fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm,999.);
      fastjet::ClusterSequence cs(particles, jet_def);
      std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets());

      if(tempJets.size() == 0)
      {
         zgs_.push_back(vector<double>());
         continue;
      }

      PseudoJet CurrentJet = tempJets[0];
      PseudoJet Part1, Part2;

      std::vector<double> z;

      while(CurrentJet.has_parents(Part1, Part2))
      {
         if(CurrentJet.pt2() <= 0)
            break;

         double DeltaR = std::sqrt(Part1.squared_distance(Part2));
         if(DeltaR < rcut_)
            break;

         double PT1 = Part1.pt();
         double PT2 = Part2.pt();
         double zg = -1;

         if(PT1 + PT2 > 0)
            zg = min(PT1, PT2) / (PT1 + PT2);
         else
            break;

         double Threshold = zcut_ * std::pow(DeltaR / r0_, beta_);

         if(zg >= Threshold)   // yay
            z.push_back(zg);

         if(PT1 > PT2)
            CurrentJet = Part1;
         else
            CurrentJet = Part2;
      }

      zgs_.push_back(z);
   }
}

std::vector<double> softDropCounter::calculateNSD(double Kappa)
{
   std::vector<double> Result;

   for(int i = 0; i < (int)zgs_.size(); i++)
   {
      double Total = 0;
      for(int j = 0; j < (int)zgs_[i].size(); j++)
         Total = Total + pow(zgs_[i][j], Kappa);
      Result.push_back(Total);
   }

   return Result;
}

#endif
