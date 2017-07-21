#ifndef treeWriter_h
#define treeWriter_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

//ROOT stuff
#include <TTree.h>

#include "fastjet/PseudoJet.hh"

using namespace std;

//---------------------------------------------------------------
// Description
// This class writes a root tree
// Only accepts vectors of the following types: int, double, fastjet::PseudoJet
// In case of PseudoJet it will store pt, eta, phi and mass as separate vectors
// in the output tree
// Authors: Y. Chen, M. Verweij
//---------------------------------------------------------------

class treeWriter {

 private :
  TTree* treeOut_;
  const char *treeName_;
  std::map<std::string,std::vector<int>  > intMaps_;
  std::map<std::string,std::vector<double>  > doubleMaps_;
      
 public :
  treeWriter(const char *treeName = "treeOut") :
    treeName_(treeName)
  {
    treeOut_ = new TTree(treeName_,"JetToyHI tree");
    
  }

  TTree *getTree() const {return treeOut_;};

  void setTreeName(const char *c) {treeName_ = c; treeOut_->SetName(c);}

  void fillTree() {treeOut_->Fill();}

  void addJetCollection(std::string name, std::vector<fastjet::PseudoJet> v) {

    //we are storing the pt, eta, phi and mass of the jets
    std::vector<double> pt;  pt.reserve(v.size());
    std::vector<double> eta; eta.reserve(v.size());
    std::vector<double> phi; phi.reserve(v.size());
    std::vector<double> m;   m.reserve(v.size());
    for( fastjet::PseudoJet jet : v ) {
      pt.push_back(jet.pt());
      eta.push_back(jet.eta());
      phi.push_back(jet.phi());
      m.push_back(jet.m());
    }
    addDoubleCollection(name+"Pt",pt);
    addDoubleCollection(name+"Eta",eta);
    addDoubleCollection(name+"Phi",phi);
    addDoubleCollection(name+"M",m);
  }

  void addDoubleCollection(std::string name, std::vector<double> v) {
    doubleMaps_[name] = v;
    bookBranchDoubleVec(name);
  }

  void addIntCollection(std::string name, std::vector<int> v) {
    intMaps_[name] = v;
    bookBranchIntVec(name);
  }

  void bookBranchDoubleVec(std::string name) {
    if(!treeOut_->GetBranch(name.c_str())) treeOut_->Branch(name.c_str(),&doubleMaps_[name]);
  }

  void bookBranchIntVec(std::string name) {
    if(!treeOut_->GetBranch(name.c_str())) treeOut_->Branch(name.c_str(),&intMaps_[name]);
  }

};
#endif
