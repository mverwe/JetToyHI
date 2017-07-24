#ifndef treeWriter_h
#define treeWriter_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

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
  std::map<std::string,std::vector<int>  > intMaps_;
  std::map<std::string,std::vector<double>  > doubleMaps_;
      
 public :
  treeWriter(const char *treeName = "treeOut");
  TTree *getTree() const;
  void setTreeName(const char *c);
  void fillTree();
  void addCollection(std::string name, const jetCollection &c);
  void addCollection(std::string name, const std::vector<fastjet::PseudoJet> &v);
  void addCollection(std::string name, const std::vector<double> &v);
  void addCollection(std::string name, const std::vector<int> &v);
  void addJetCollection(std::string name, const jetCollection &c);
  void addJetCollection(std::string name, const std::vector<fastjet::PseudoJet> v);
  void addDoubleCollection(std::string name, const std::vector<double> v);
  void addIntCollection(std::string name, const std::vector<int> v);
  void bookBranchDoubleVec(std::string name);
  void bookBranchIntVec(std::string name);
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

void treeWriter::addCollection(std::string name, const jetCollection &c)
{
   addJetCollection(name, c);
}

void treeWriter::addCollection(std::string name, const std::vector<fastjet::PseudoJet> &v)
{
   addJetCollection(name, v);
}

void treeWriter::addCollection(std::string name, const std::vector<double> &v)
{
   addDoubleCollection(name, v);
}

void treeWriter::addCollection(std::string name, const std::vector<int> &v)
{
   addIntCollection(name, v);
}

void treeWriter::addJetCollection(std::string name, const jetCollection &c)
{
   addJetCollection(name, c.getJet());

   std::vector<std::string> doubleKeys = c.getListOfKeysDouble();
   for(std::string tag: doubleKeys)
      addDoubleCollection(tag, c.getVectorDouble(tag));

   std::vector<std::string> intKeys = c.getListOfKeysInt();
   for(std::string tag: intKeys)
      addIntCollection(tag, c.getVectorInt(tag));
}

void treeWriter::addJetCollection(std::string name, const std::vector<fastjet::PseudoJet> v)
{
   //we are storing the pt, eta, phi and mass of the jets
   std::vector<double> pt;    pt.reserve(v.size());
   std::vector<double> eta;   eta.reserve(v.size());
   std::vector<double> phi;   phi.reserve(v.size());
   std::vector<double> m;     m.reserve(v.size());
   std::vector<double> area;    area.reserve(v.size());

   for(const fastjet::PseudoJet jet: v)
   {
      pt.push_back(jet.pt());
      eta.push_back(jet.eta());
      phi.push_back(jet.phi());
      m.push_back(jet.m());
      if(jet.has_area()) area.push_back(jet.area());
      else area.push_back(-1.);
   }

   addDoubleCollection(name + "Pt",  pt);
   addDoubleCollection(name + "Eta", eta);
   addDoubleCollection(name + "Phi", phi);
   addDoubleCollection(name + "M",   m);
   addDoubleCollection(name + "Area",   area);
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



#endif
