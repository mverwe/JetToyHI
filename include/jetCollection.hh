#ifndef jetCollection_HH
#define jetCollection_HH

#include <iostream>
#include <vector>

#include "fastjet/PseudoJet.hh"

using namespace std;

//---------------------------------------------------------------
// Description
// This class contains jet momentum and associated quantities for the whole event
// Author: Y. Chen
//---------------------------------------------------------------

class jetCollection
{
private:
   std::vector<fastjet::PseudoJet> p_;
   std::map<std::string, std::vector<double>> doublemap_;
   std::map<std::string, std::vector<int>> intmap_;
public:
   jetCollection(const std::vector<fastjet::PseudoJet> &p);
   ~jetCollection();
   void setJet(const std::vector<fastjet::PseudoJet> &v);
   std::vector<fastjet::PseudoJet> getJet() const;
   std::vector<double> getVectorDouble(string tag) const;
   std::vector<int> getVectorInt(string tag) const;
   void addVector(string tag, std::vector<double> v);
   void addVector(string tag, std::vector<int> v);
   std::vector<std::string> getListOfKeysDouble() const;
   std::vector<std::string> getListOfKeysInt() const;
};

jetCollection::jetCollection(const std::vector<fastjet::PseudoJet> &p)
   : p_(p)
{
}

jetCollection::~jetCollection()
{
}
   
void jetCollection::setJet(const std::vector<fastjet::PseudoJet> &v)
{
   p_ = v;
}
   
std::vector<fastjet::PseudoJet> jetCollection::getJet() const
{
   return p_;
}
   
std::vector<double> jetCollection::getVectorDouble(string tag) const
{
   if(doublemap_.find(tag) == doublemap_.end())
      return std::vector<double>();
   return doublemap_.at(tag);
}

std::vector<int> jetCollection::getVectorInt(string tag) const
{
   if(intmap_.find(tag) == intmap_.end())
      return std::vector<int>();
   return intmap_.at(tag);
}

void jetCollection::addVector(string tag, std::vector<double> v)
{
   if(intmap_.find(tag) != intmap_.end())
      intmap_.erase(tag);
   doublemap_[tag] = v;
}

void jetCollection::addVector(string tag, std::vector<int> v)
{
   if(doublemap_.find(tag) != doublemap_.end())
      doublemap_.erase(tag);
   intmap_[tag] = v;
}

std::vector<std::string> jetCollection::getListOfKeysDouble() const
{
   std::vector<std::string> Result;
   for(std::map<std::string, std::vector<double>>::const_iterator iter = doublemap_.begin(); iter != doublemap_.end(); iter++)
      Result.push_back(iter->first);
   return Result;
}

std::vector<std::string> jetCollection::getListOfKeysInt() const
{
   std::vector<std::string> Result;
   for(std::map<std::string, std::vector<int>>::const_iterator iter = intmap_.begin(); iter != intmap_.end(); iter++)
      Result.push_back(iter->first);
   return Result;
}

#endif
