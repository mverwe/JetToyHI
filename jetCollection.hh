#ifndef jetCollection_HH
#define jetCollection_HH

#include <iostream>
#include <vector>

#include "fastjet/PseudoJet.hh"

class jetCollection
{
public:
   std::vector<fastjet::PseudoJet> p_;
   std::map<std::string, std::vector<double>> doublemap_;
   std::map<std::string, std::vector<int>> intmap_;
public:
   jetCollection(const std::vector<fastjet::PseudoJet> &p);
   ~jetCollection();
   std::vector<fastjet::PseudoJet> getJet();
   std::vector<double> getVectorDouble(string tag);
   std::vector<int> getVectorInt(string tag);
   void addVector(string tag, std::vector<double> v);
   void addVector(string tag, std::vector<int> v);
   std::vector<std::string> getListOfKeysDouble();
   std::vector<std::string> getListOfKeysInt();
};

jetCollection::jetCollection(const std::vector<fastjet::PseudoJet> &p)
   : p_(p)
{
}

jetCollection::~jetCollection()
{
}
   
std::vector<fastjet::PseudoJet> jetCollection::getJet()
{
   return p_;
}
   
std::vector<double> jetCollection::getVectorDouble(string tag)
{
   if(doublemap_.find(tag) == doublemap_.end())
      return std::vector<double>();
   return doublemap_[tag];
}

std::vector<int> jetCollection::getVectorInt(string tag)
{
   if(intmap_.find(tag) == intmap_.end())
      return std::vector<int>();
   return intmap_[tag];
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

std::vector<std::string> jetCollection::getListOfKeysDouble()
{
   std::vector<std::string> Result;
   for(std::map<std::string, std::vector<double>>::iterator iter = doublemap_.begin(); iter != doublemap_.end(); iter++)
      Result.push_back(iter->first);
   return Result;
}

std::vector<std::string> jetCollection::getListOfKeysInt()
{
   std::vector<std::string> Result;
   for(std::map<std::string, std::vector<int>>::iterator iter = intmap_.begin(); iter != intmap_.end(); iter++)
      Result.push_back(iter->first);
   return Result;
}

#endif
