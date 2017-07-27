#ifndef skSubtractor_h
#define skSubtractor_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "fastjet/contrib/SoftKiller.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs soft killer on the full event
// Author: M. Verweij
//---------------------------------------------------------------

class skSubtractor {

private :
  double gridSize_;
  double rapMax_;
  double ptThreshold_;
  std::vector<fastjet::PseudoJet> fjInputs_;

  contrib::SoftKiller subtractor_;

  
public :
  skSubtractor(double gridSize = 0.4, double rapMax = 3.) :
    gridSize_(gridSize),
    rapMax_(rapMax),
    ptThreshold_()
  {
    //init SoftKiller subtractor
    subtractor_ = contrib::SoftKiller(rapMax_,gridSize_);
  }

  void setInputParticles(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }

  double getPtThreshold()  const { return ptThreshold_; }
  
  std::vector<fastjet::PseudoJet> doSubtraction() {

    std::vector<fastjet::PseudoJet> skEvent;
    subtractor_.apply(fjInputs_, skEvent, ptThreshold_);

    // std::cout << "# Ran the following soft killer: " << subtractor_.description() << std::endl;
    
    return skEvent;
  }
};

#endif
