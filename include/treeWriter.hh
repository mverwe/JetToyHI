#ifndef treeWriter_h
#define treeWriter_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <map>

#include "TTree.h"

#include "fastjet/PseudoJet.hh"

#include "jetCollection.hh"

//---------------------------------------------------------------
// Description
// This class writes a root tree
// Only accepts vectors of the following types: int, double, fastjet::PseudoJet
// In case of PseudoJet it will store pt, eta, phi and mass as separate vectors
// in the output tree
// Authors: Y. Chen, M. Verweij
//---------------------------------------------------------------

class treeWriter
{
private :
  TTree* treeOut_;
  const char *treeName_;
  std::map<std::string,std::vector<bool>  > boolMaps_;
  std::map<std::string,std::vector<int>  > intMaps_;
  std::map<std::string,std::vector<double>  > doubleMaps_;
  std::map<std::string,std::vector<std::vector<double>>> doubleVectorMaps_;
  std::map<std::string,std::vector<std::vector<int>>> intVectorMaps_;

public :
  treeWriter(const char *treeName = "treeOut");
  TTree *getTree() const;
  void setTreeName(const char *c);
  void fillTree();
  void addCollection(std::string name, const jetCollection &c, bool writeConst = false);
  void addCollection(std::string name, const std::vector<fastjet::PseudoJet> &v, bool writeConst = false);
  void addCollection(std::string name, const std::vector<double> &v);
  void addCollection(std::string name, const std::vector<int> &v);
  void addCollection(std::string name, const std::vector<bool> &v);
  void addJetCollection(std::string name, const jetCollection &c, bool writeConst = false);
  void addJetCollection(std::string name, const std::vector<fastjet::PseudoJet> v, bool writeConst = false);

  void addPartonCollection(std::string name, const std::vector<fastjet::PseudoJet> v);
  
  void addDoubleCollection(std::string name, const std::vector<double> v);
  void addIntCollection(std::string name, const std::vector<int> v);
  void addBoolCollection(std::string name, const std::vector<bool> v);
  void bookBranchDoubleVec(std::string name);
  void bookBranchIntVec(std::string name);
  void bookBranchBoolVec(std::string name);

  void addDoubleVectorCollection(std::string name, const std::vector<std::vector<double>> v);
  void addIntVectorCollection(std::string name, const std::vector<std::vector<int>> v);
  void bookBranchDoubleVectorVec(std::string name);
  void bookBranchIntVectorVec(std::string name);

};

treeWriter::treeWriter(const char *treeName)
  : treeName_(treeName)
{
  treeOut_ = new TTree(treeName_,"JetToyHI tree");
}

TTree *treeWriter::getTree() const
{
  return treeOut_;
}

void treeWriter::setTreeName(const char *c)
{
  treeName_ = c;
  treeOut_->SetName(c);
}

void treeWriter::fillTree()
{
  treeOut_->Fill();
}

void treeWriter::addCollection(std::string name, const jetCollection &c, bool writeConst)
{
  addJetCollection(name, c, writeConst);
}

void treeWriter::addCollection(std::string name, const std::vector<fastjet::PseudoJet> &v, bool writeConst)
{
  addJetCollection(name, v, writeConst);
}

void treeWriter::addCollection(std::string name, const std::vector<double> &v)
{
  addDoubleCollection(name, v);
}

void treeWriter::addCollection(std::string name, const std::vector<int> &v)
{
  addIntCollection(name, v);
}

void treeWriter::addCollection(std::string name, const std::vector<bool> &v)
{
  addBoolCollection(name, v);
}


void treeWriter::addJetCollection(std::string name, const jetCollection &c, bool writeConst)
{
  addJetCollection(name, c.getJet(), writeConst);

  std::vector<std::string> jetKeys = c.getListOfKeysJet();
  for(std::string tag: jetKeys)
    addJetCollection(tag, c.getVectorJet(tag));

  std::vector<std::string> doubleKeys = c.getListOfKeysDouble();
  for(std::string tag: doubleKeys)
    addDoubleCollection(tag, c.getVectorDouble(tag));

  std::vector<std::string> doubledoubleKeys = c.getListOfKeysDoubleDouble();
  for(std::string tag: doubledoubleKeys)
     addDoubleVectorCollection(tag, c.getVectorDoubleDouble(tag));

  std::vector<std::string> intintKeys = c.getListOfKeysIntInt();
  for(std::string tag: intintKeys)
     addIntVectorCollection(tag, c.getVectorIntInt(tag));

  std::vector<std::string> intKeys = c.getListOfKeysInt();
  for(std::string tag: intKeys)
    addIntCollection(tag, c.getVectorInt(tag));
}

void treeWriter::addJetCollection(std::string name, const std::vector<fastjet::PseudoJet> v, bool writeConst)
{
  //we are storing the pt, eta, phi and mass of the jets
  std::vector<double> pt;    pt.reserve(v.size());
  std::vector<double> eta;   eta.reserve(v.size());
  std::vector<double> phi;   phi.reserve(v.size());
  std::vector<double> m;     m.reserve(v.size());
  std::vector<double> area;  area.reserve(v.size());

  std::vector<std::vector<double>> constPt;
  std::vector<std::vector<double>> constEta;
  std::vector<std::vector<double>> constPhi;
  std::vector<std::vector<double>> constM;
  //if(writeConst) {
  //  constPt.clear();
  //  constPt.reserve(v.size());
  //}
  for(const fastjet::PseudoJet jet: v) {
    pt.push_back(jet.pt());
    eta.push_back(jet.eta());
    phi.push_back(jet.phi());
    m.push_back(jet.m());
    if(jet.has_area()) area.push_back(jet.area());
    else area.push_back(-1.);

    if(writeConst && jet.has_valid_cluster_sequence()) {
      //get constituents of jet
      std::vector<double> ptConst;
      std::vector<double> etaConst;
      std::vector<double> phiConst;
      std::vector<double> mConst;
      
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
      for(const fastjet::PseudoJet& p : particles) {
        ptConst.push_back(p.pt());
        etaConst.push_back(p.eta());
        phiConst.push_back(p.phi());
        mConst.push_back(p.m());
      }
      constPt.push_back(ptConst);
      constEta.push_back(etaConst);
      constPhi.push_back(phiConst);
      constM.push_back(mConst);
    }
  }

  addDoubleCollection(name + "Pt",  pt);
  addDoubleCollection(name + "Eta", eta);
  addDoubleCollection(name + "Phi", phi);
  addDoubleCollection(name + "M",   m);
  addDoubleCollection(name + "Area",   area);

  if(writeConst) {
    //addDoubleVectorCollection(name + "ConstPt", constPt);
    std::string branchName = name + "ConstPt";
    doubleVectorMaps_[branchName] = constPt;
    if(!treeOut_->GetBranch(branchName.c_str()))
      treeOut_->Branch(branchName.c_str(),&doubleVectorMaps_[branchName]);
    
    branchName = name + "ConstEta";
    doubleVectorMaps_[branchName] = constEta;
    if(!treeOut_->GetBranch(branchName.c_str()))
      treeOut_->Branch(branchName.c_str(),&doubleVectorMaps_[branchName]);
    
    branchName = name + "ConstPhi";
    doubleVectorMaps_[branchName] = constPhi;
    if(!treeOut_->GetBranch(branchName.c_str()))
      treeOut_->Branch(branchName.c_str(),&doubleVectorMaps_[branchName]);

    branchName = name + "ConstM";
    doubleVectorMaps_[branchName] = constM;
    if(!treeOut_->GetBranch(branchName.c_str()))
      treeOut_->Branch(branchName.c_str(),&doubleVectorMaps_[branchName]);
    
  }
}

void treeWriter::addPartonCollection(std::string name, const std::vector<fastjet::PseudoJet> v)
{
  //we are storing the pt, eta, phi, mass and pdg of partons
  std::vector<double> pt;    pt.reserve(v.size());
  std::vector<double> eta;   eta.reserve(v.size());
  std::vector<double> phi;   phi.reserve(v.size());
  std::vector<double> m;     m.reserve(v.size());
  std::vector<int>    pdg;   pdg.reserve(v.size());

  for(const fastjet::PseudoJet p: v) {
    pt.push_back(p.pt());
    eta.push_back(p.eta());
    phi.push_back(p.phi());
    m.push_back(p.m());
    const int & pdgid = p.user_info<PU14>().pdg_id();
    pdg.push_back(pdgid);
  }

  addDoubleCollection(name + "Pt",  pt);
  addDoubleCollection(name + "Eta", eta);
  addDoubleCollection(name + "Phi", phi);
  addDoubleCollection(name + "M",   m);
  addIntCollection(name + "PDG",   pdg);

}

void treeWriter::addDoubleCollection(std::string name, const std::vector<double> v)
{
  doubleMaps_[name] = v;
  bookBranchDoubleVec(name);
}

void treeWriter::addIntCollection(std::string name, const std::vector<int> v)
{
  intMaps_[name] = v;
  bookBranchIntVec(name);
}

void treeWriter::addBoolCollection(std::string name, const std::vector<bool> v)
{
  boolMaps_[name] = v;
  bookBranchBoolVec(name);
}

void treeWriter::bookBranchDoubleVec(std::string name)
{
  if(!treeOut_->GetBranch(name.c_str()))
    treeOut_->Branch(name.c_str(),&doubleMaps_[name]);
}

void treeWriter::bookBranchIntVec(std::string name)
{
  if(!treeOut_->GetBranch(name.c_str()))
    treeOut_->Branch(name.c_str(),&intMaps_[name]);
}

void treeWriter::bookBranchBoolVec(std::string name)
{
  if(!treeOut_->GetBranch(name.c_str()))
    treeOut_->Branch(name.c_str(),&boolMaps_[name]);
}


void treeWriter::addDoubleVectorCollection(std::string name, const std::vector<std::vector<double> > v)
{
  doubleVectorMaps_[name] = v;
  bookBranchDoubleVectorVec(name);
}

void treeWriter::addIntVectorCollection(std::string name, const std::vector<std::vector<int> > v)
{
  intVectorMaps_[name] = v;
  bookBranchIntVectorVec(name);
}


void treeWriter::bookBranchDoubleVectorVec(std::string name)
{
  if(!treeOut_->GetBranch(name.c_str()))
    treeOut_->Branch(name.c_str(), &doubleVectorMaps_[name]);
}

void treeWriter::bookBranchIntVectorVec(std::string name)
{
  if(!treeOut_->GetBranch(name.c_str()))
    treeOut_->Branch(name.c_str(), &intVectorMaps_[name]);
}


#endif
