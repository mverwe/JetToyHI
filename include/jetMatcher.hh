#ifndef jetMatcher_h
#define jetMatcher_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"

#include "jetCollection.hh"

//---------------------------------------------------------------
// Description
// This class matches two collections of jets using bijective algorithm
// Author: M. Verweij, Y. Chen
//---------------------------------------------------------------

class jetMatcher {

 private :
  std::vector<fastjet::PseudoJet> fjBase_;   //base jets
  std::vector<fastjet::PseudoJet> fjTag_;    //jets to match to base

  std::vector<int>                fjBaseMatchIds_; //base jet -> tag jet ID
  std::vector<int>                fjTagMatchIds_;  //tag jet -> base jet ID

  double                          maxDist_;  //max distance for matching (to limit CPU time)

 public :
  jetMatcher(double maxDist = 0.4);
  void setBaseJets(std::vector<fastjet::PseudoJet> v);
  void setTagJets(std::vector<fastjet::PseudoJet> v);
  void setBaseJets(jetCollection c);
  void setTagJets(jetCollection c);
  void setMaxDist(double d);

  std::vector<int> getBaseMatchIds() const;
  std::vector<int> getTagMatchIds()  const;

  void matchJets();
  std::vector<fastjet::PseudoJet> getTagJetsOrderedToBase();
  std::vector<fastjet::PseudoJet> getBaseJetsOrderedToTag();
  std::vector<fastjet::PseudoJet> reorderedToBase(std::vector<fastjet::PseudoJet> v);
  std::vector<fastjet::PseudoJet> reorderedToTag(std::vector<fastjet::PseudoJet> v);
  void reorderedToBase(jetCollection &c);
  void reorderedToTag(jetCollection &c);
  std::vector<double> reorderedToBase(std::vector<double> v);
  std::vector<double> reorderedToTag(std::vector<double> v);
  std::vector<int> reorderedToBase(std::vector<int> v);
  std::vector<int> reorderedToTag(std::vector<int> v);
};

jetMatcher::jetMatcher(double maxDist)
   : maxDist_(maxDist)
{
}

void jetMatcher::setBaseJets(std::vector<fastjet::PseudoJet> v)
{
   fjBase_  = v;
}

void jetMatcher::setTagJets(std::vector<fastjet::PseudoJet> v)
{
   fjTag_   = v;
}

void jetMatcher::setBaseJets(jetCollection c)
{
   fjBase_  = c.getJet();
}

void jetMatcher::setTagJets(jetCollection c)
{
   fjTag_   = c.getJet();
}

void jetMatcher::setMaxDist(double d)
{
   maxDist_ = d;
}

std::vector<int> jetMatcher::getBaseMatchIds() const
{
   return fjBaseMatchIds_;
}

std::vector<int> jetMatcher::getTagMatchIds()  const
{
   return fjTagMatchIds_;
}

void jetMatcher::matchJets()
{
   int nJets1 = (int)fjBase_.size();
   int nJets2 = (int)fjTag_.size();

   fjBaseMatchIds_.resize(fjBase_.size());
   fjTagMatchIds_.resize(fjTag_.size());
   std::fill(fjBaseMatchIds_.begin(),fjBaseMatchIds_.end(),-1);
   std::fill(fjTagMatchIds_.begin(),fjTagMatchIds_.end(),-1);

   std::vector<int> matchIndex1;
   matchIndex1.resize(nJets2+1);
   std::fill(matchIndex1.begin(), matchIndex1.end(), -1);

   std::vector<int> matchIndex2;
   matchIndex2.resize(nJets1+1);
   std::fill(matchIndex2.begin(), matchIndex2.end(), -1);

   std::vector<int> iFlag;
   iFlag.resize(nJets1*nJets2);
   std::fill(iFlag.begin(), iFlag.end(), 0);

   //Find the clostest match jet to base jet
   Double_t dist = maxDist_;
   for (int i = 0; i < nJets1; i++) {
      fastjet::PseudoJet jet1 = fjBase_[i];
      if(fabs(jet1.pt())<1e-6) continue; //remove ghosts

      dist = maxDist_;
      for (int j = 0; j < nJets2; j++) {
         fastjet::PseudoJet jet2 = fjTag_[j];
         if(fabs(jet2.pt())<1e-6) continue; //remove ghosts

         double dR = jet1.delta_R(jet2);
         if(dR<dist && dR<maxDist_){
            matchIndex2[i]=j;
            dist = dR;
         }
      }//jet2 loop
      if(matchIndex2[i]>=0) {
         iFlag[i*nJets2+matchIndex2[i]]+=1;//j closest to i
      }
   }//jet1 loop

   //other way around
   for (int j = 0; j < nJets2; j++) {
      fastjet::PseudoJet jet2 = fjTag_[j];
      if(fabs(jet2.pt())<1e-6) continue; //remove ghosts

      dist = maxDist_;
      for (int i = 0; i < nJets1; i++) {
         fastjet::PseudoJet jet1 = fjBase_[i];
         if(fabs(jet1.pt())<1e-6) continue; //remove ghosts

         double dR = jet1.delta_R(jet2);
         if(dR<dist && dR<maxDist_){
            matchIndex1[j]=i;
            dist = dR;
         }
      }
      if(matchIndex1[j]>=0) {
         iFlag[matchIndex1[j]*nJets2+j]+=2;//i closest to j
      }
   }//jet2 loop

   // check for "true" correlations
   for (int i = 0; i < nJets1; i++) {
      fastjet::PseudoJet jet1 = fjBase_[i];
      if(fabs(jet1.pt())<1e-6) continue; //remove ghosts

      dist = maxDist_;
      for (int j = 0; j < nJets2; j++) {
         fastjet::PseudoJet jet2 = fjTag_[j];
         if(fabs(jet2.pt())<1e-6) continue; //remove ghosts

         // do we have a unique correlation?
         if(iFlag[i*nJets2+j]==3) {
            fjBaseMatchIds_[i] = j;
            fjTagMatchIds_[j] = i;
         }
      }
   }

}

std::vector<fastjet::PseudoJet> jetMatcher::getTagJetsOrderedToBase()
{
   std::vector<fastjet::PseudoJet> tagReordered;
   tagReordered.resize(fjBase_.size());
   for (unsigned int j = 0; j < fjTag_.size(); j++) {
      if(fjTagMatchIds_[j]<0) continue;
      tagReordered[fjTagMatchIds_[j]] = fjTag_[j];
   }
   return tagReordered;
}

std::vector<fastjet::PseudoJet> jetMatcher::getBaseJetsOrderedToTag()
{
   std::vector<fastjet::PseudoJet> baseReordered;
   baseReordered.resize(fjTag_.size());
   for (unsigned int i = 0; i < fjBase_.size(); i++) {
      if(fjBaseMatchIds_[i]<0) continue;
      baseReordered[fjBaseMatchIds_[i]] = fjBase_[i];
   }
   return baseReordered;
}

std::vector<fastjet::PseudoJet> jetMatcher::reorderedToBase(std::vector<fastjet::PseudoJet> v)
{
   std::vector<fastjet::PseudoJet> vecReordered;
   vecReordered.resize(fjBase_.size());
   if(v.size() != fjTag_.size()) {
      std::cout << "WARNING: vector size not compatible. Not ordering" << std::endl;
      return vecReordered;
   }
   for (unsigned int i = 0; i < v.size(); i++) {
      if(fjTagMatchIds_[i]<0) continue;
      vecReordered[fjTagMatchIds_[i]] = v[i];
   }
   return vecReordered;
}

std::vector<fastjet::PseudoJet> jetMatcher::reorderedToTag(std::vector<fastjet::PseudoJet> v)
{
   std::vector<fastjet::PseudoJet> vecReordered;
   vecReordered.resize(fjTag_.size());
   if(v.size() != fjBase_.size()) {
      std::cout << "WARNING: vector size not compatible. Not ordering" << std::endl;
      return vecReordered;
   }
   for (unsigned int i = 0; i < v.size(); i++) {
      if(fjBaseMatchIds_[i]<0) continue;
      vecReordered[fjBaseMatchIds_[i]] = v[i];
   }
   return vecReordered;
}

void jetMatcher::reorderedToBase(jetCollection &c)
{
   c.setJet(reorderedToBase(c.getJet()));

   std::vector<std::string> doubleKeys = c.getListOfKeysDouble();
   for(std::string tag: doubleKeys)
      c.addVector(tag, reorderedToBase(c.getVectorDouble(tag)));

   std::vector<std::string> intKeys = c.getListOfKeysInt();
   for(std::string tag: intKeys)
      c.addVector(tag, reorderedToBase(c.getVectorInt(tag)));
};

void jetMatcher::reorderedToTag(jetCollection &c)
{
   c.setJet(reorderedToTag(c.getJet()));

   std::vector<std::string> doubleKeys = c.getListOfKeysDouble();
   for(std::string tag: doubleKeys)
      c.addVector(tag, reorderedToTag(c.getVectorDouble(tag)));

   std::vector<std::string> intKeys = c.getListOfKeysInt();
   for(std::string tag: intKeys)
      c.addVector(tag, reorderedToTag(c.getVectorInt(tag)));
};

std::vector<double> jetMatcher::reorderedToBase(std::vector<double> v)
{
   std::vector<double> vecReordered;
   vecReordered.resize(fjBase_.size());
   if(v.size() != fjTag_.size()) {
      std::cout << "WARNING: vector size not compatible. Not ordering" << std::endl;
      return vecReordered;
   }
   for (unsigned int i = 0; i < v.size(); i++) {
      if(fjTagMatchIds_[i]<0) continue;
      vecReordered[fjTagMatchIds_[i]] = v[i];
   }
   return vecReordered;
}

std::vector<double> jetMatcher::reorderedToTag(std::vector<double> v)
{
   std::vector<double> vecReordered;
   vecReordered.resize(fjTag_.size());
   if(v.size() != fjBase_.size()) {
      std::cout << "WARNING: vector size not compatible. Not ordering" << std::endl;
      return vecReordered;
   }
   for (unsigned int i = 0; i < v.size(); i++) {
      if(fjBaseMatchIds_[i]<0) continue;
      vecReordered[fjBaseMatchIds_[i]] = v[i];
   }
   return vecReordered;
}

std::vector<int> jetMatcher::reorderedToBase(std::vector<int> v)
{
   std::vector<int> vecReordered;
   vecReordered.resize(fjBase_.size());
   if(v.size() != fjTag_.size()) {
      std::cout << "WARNING: vector size not compatible. Not ordering" << std::endl;
      return vecReordered;
   }
   for (unsigned int i = 0; i < v.size(); i++) {
      if(fjTagMatchIds_[i]<0) continue;
      vecReordered[fjTagMatchIds_[i]] = v[i];
   }
   return vecReordered;
}

std::vector<int> jetMatcher::reorderedToTag(std::vector<int> v)
{
   std::vector<int> vecReordered;
   vecReordered.resize(fjTag_.size());
   if(v.size() != fjBase_.size()) {
      std::cout << "WARNING: vector size not compatible. Not ordering" << std::endl;
      return vecReordered;
   }
   for (unsigned int i = 0; i < v.size(); i++) {
      if(fjBaseMatchIds_[i]<0) continue;
      vecReordered[fjBaseMatchIds_[i]] = v[i];
   }
   return vecReordered;
}

#endif
