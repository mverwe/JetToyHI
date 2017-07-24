/// This file is provides a collection of helpers:
///
/// MasslessTransformer   a FastJet Transformer that sets particle masses to 0
/// Width                 Width/Girth/Angularity(1) jet shape
///

#ifndef __PU14_HELPERS_HH__
#define __PU14_HELPERS_HH__

#include "fastjet/tools/Transformer.hh"

/// \class MasslessTransformer
/// a FastJet transformer that sets the mass of the particles/jets to
/// 0 (preserving all the other info, e.g. rapidity)
class MasslessTransformer : public fastjet::Transformer{
public:
  MasslessTransformer(){};
  virtual std::string description() const{ return "makes a particle massless, conserving rapidity";}
  virtual fastjet::PseudoJet result(const fastjet::PseudoJet &jet) const{
    fastjet::PseudoJet j = jet;
    j.reset_momentum(fastjet::PtYPhiM(jet.pt(), jet.rap(), jet.phi()));
    return j;
  }
};


/// \class Width
class Width : public fastjet::FunctionOfPseudoJet<double>{
public:
  Width(){}
  virtual std::string description() const{ return "width/girth/angularity(1)";}
  virtual double result(const fastjet::PseudoJet &jet) const;
};
  

#endif  // __PU14_HELPERS_HH__
