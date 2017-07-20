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

class treeWriter {

 private :
  TTree* treeOut_;
  const char *treeName_;
  std::map<std::string,std::vector<int>  > intMaps_;
  std::map<std::string,std::vector<double>  > doubleMaps_;
  
  
  std::vector<std::map<std::string,std::vector<fastjet::PseudoJet> > > jetMaps_;
    
 public :
  treeWriter() {
    treeOut_ = new TTree("treeOut","JetToyHI tree");
    
  }

  TTree *getTree() const {return treeOut_;};

  void setTreeName(const char *c) {treeName_ = c; treeOut_->SetName(c);}

  void fillTree() {treeOut_->Fill();}

  void addJetCollection(std::string name, std::vector<fastjet::PseudoJet> v) {
    std::map<std::string,std::vector<fastjet::PseudoJet> > myMap;
    for(fastjet::PseudoJet val : v)
      myMap[name].push_back(val);
    jetMaps_.push_back(myMap);
  }

  void addDoubleCollection(std::string name, std::vector<double> v) {
    doubleMaps_[name] = v;
    bookBranchVec(name);
  }

  void addIntCollection(std::string name, std::vector<int> v) {
    intMaps_[name] = v;
  }

  void bookBranchVec(std::string name) { //, std::vector<double> var) {
    if(!treeOut_->GetBranch(name.c_str())) treeOut_->Branch(name.c_str(),&doubleMaps_[name]);
  }

};
#endif
