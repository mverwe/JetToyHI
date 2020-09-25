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
   std::map<std::string, std::vector<fastjet::PseudoJet>> jetmap_;
   std::map<std::string, std::vector<double>> doublemap_;
   std::map<std::string, std::vector<int>> intmap_;
   std::map<std::string, std::vector<std::vector<double>>> doubledoublemap_;
   std::map<std::string, std::vector<std::vector<int>>> intintmap_;
public:
   jetCollection(const std::vector<fastjet::PseudoJet> &p);
   ~jetCollection();
   void setJet(const std::vector<fastjet::PseudoJet> &v);
   std::vector<fastjet::PseudoJet> getJet() const;
   std::vector<fastjet::PseudoJet> getVectorJet(string tag) const;
   std::vector<double> getVectorDouble(string tag) const;
   std::vector<std::vector<double>> getVectorDoubleDouble(string tag) const;
   std::vector<int> getVectorInt(string tag) const;
   std::vector<std::vector<int>> getVectorIntInt(string tag) const;
   void addVector(string tag, std::vector<fastjet::PseudoJet> v);
   void addVector(string tag, std::vector<double> v);
   void addVector(string tag, std::vector<std::vector<double>> v);
   void addVector(string tag, std::vector<int> v);
   void addVector(string tag, std::vector<std::vector<int>> v);
   std::vector<std::string> getListOfKeysJet() const;
   std::vector<std::string> getListOfKeysDouble() const;
   std::vector<std::string> getListOfKeysDoubleDouble() const;
   std::vector<std::string> getListOfKeysInt() const;
   std::vector<std::string> getListOfKeysIntInt() const;
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
   
std::vector<fastjet::PseudoJet> jetCollection::getVectorJet(string tag) const
{
   if(jetmap_.find(tag) == jetmap_.end())
      return std::vector<fastjet::PseudoJet>();
   return jetmap_.at(tag);
}

std::vector<double> jetCollection::getVectorDouble(string tag) const
{
   if(doublemap_.find(tag) == doublemap_.end())
      return std::vector<double>();
   return doublemap_.at(tag);
}

std::vector<std::vector<double>> jetCollection::getVectorDoubleDouble(string tag) const
{
   if(doubledoublemap_.find(tag) == doubledoublemap_.end())
      return std::vector<std::vector<double>>();
   return doubledoublemap_.at(tag);
}

std::vector<int> jetCollection::getVectorInt(string tag) const
{
   if(intmap_.find(tag) == intmap_.end())
      return std::vector<int>();
   return intmap_.at(tag);
}

std::vector<std::vector<int>> jetCollection::getVectorIntInt(string tag) const
{
   if(intintmap_.find(tag) == intintmap_.end())
      return std::vector<std::vector<int>>();
   return intintmap_.at(tag);
}


void jetCollection::addVector(string tag, std::vector<fastjet::PseudoJet> v)
{
   if(jetmap_.find(tag) != jetmap_.end())
      jetmap_.erase(tag);
   jetmap_[tag] = v;
}

void jetCollection::addVector(string tag, std::vector<int> v)
{
   if(intmap_.find(tag) != intmap_.end())
      intmap_.erase(tag);
   intmap_[tag] = v;
}

void jetCollection::addVector(string tag, std::vector<double> v)
{
   if(doublemap_.find(tag) != doublemap_.end())
      doublemap_.erase(tag);
   doublemap_[tag] = v;
}

void jetCollection::addVector(string tag, std::vector<std::vector<double>> v)
{
   if(doubledoublemap_.find(tag) != doubledoublemap_.end())
      doubledoublemap_.erase(tag);
   doubledoublemap_[tag] = v;
}

void jetCollection::addVector(string tag, std::vector<std::vector<int>> v)
{
   if(intintmap_.find(tag) != intintmap_.end())
      intintmap_.erase(tag);
   intintmap_[tag] = v;
}


std::vector<std::string> jetCollection::getListOfKeysJet() const
{
   std::vector<std::string> Result;
   for(auto iter = jetmap_.begin(); iter != jetmap_.end(); iter++)
      Result.push_back(iter->first);
   return Result;
}

std::vector<std::string> jetCollection::getListOfKeysDouble() const
{
   std::vector<std::string> Result;
   for(auto iter = doublemap_.begin(); iter != doublemap_.end(); iter++)
      Result.push_back(iter->first);
   return Result;
}

std::vector<std::string> jetCollection::getListOfKeysDoubleDouble() const
{
   std::vector<std::string> Result;
   for(auto iter = doubledoublemap_.begin(); iter != doubledoublemap_.end(); iter++)
      Result.push_back(iter->first);
   return Result;
}

std::vector<std::string> jetCollection::getListOfKeysInt() const
{
   std::vector<std::string> Result;
   for(auto iter = intmap_.begin(); iter != intmap_.end(); iter++)
      Result.push_back(iter->first);
   return Result;
}

std::vector<std::string> jetCollection::getListOfKeysIntInt() const
{
   std::vector<std::string> Result;
   for(auto iter = intintmap_.begin(); iter != intintmap_.end(); iter++)
      Result.push_back(iter->first);
   return Result;
}


#endif
