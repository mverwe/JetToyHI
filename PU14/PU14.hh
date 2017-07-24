#ifndef __PU14_HH__
#define __PU14_HH__

#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "HepPID/ParticleIDMethods.hh"

//----------------------------------------------------------------------
/// \class PU14
/// 
/// class that stores basic user information for the 2014 pileup
/// workshop.
///
class PU14 : public fastjet::PseudoJet::UserInfoBase {
public:
  /// construct the user info for the PU14 studies
  PU14(int pdg_id, int barcode = 0, int vertex = 0) {
    _pdg_id = pdg_id;
    _three_charge = HepPID::threeCharge(pdg_id);
    _barcode = barcode;
    _vertex = vertex;
  }
  
  /// returns the pdg ID of the particle
  int pdg_id() const {return _pdg_id;}
  
  /// returns its charge
  double charge() const {return _three_charge / 3.0;}

  /// returns three times its charge (an integer)
  int three_charge() const {return _three_charge;}

  /// returns the vertex number (0 is primary hard vertex)
  int vertex_number() const {return _vertex;}

  /// returns the vertex index (0 is primary hard vertex)
  int vertex() const {return _vertex;}

private:
  int _pdg_id, _three_charge;
  int _barcode;
  int _vertex;
};


std::ostream & operator<<(std::ostream &o , const fastjet::PseudoJet & p);


/// returns a Selector that is true for charged particles
fastjet::Selector SelectorIsCharged();
/// alternative name for charged-particle selector
inline 
fastjet::Selector SelectorCharged() {return SelectorIsCharged();}

/// returns a Selector that is true for particles that come from the
/// specified vertex
fastjet::Selector SelectorVertexNumber(int i);

/// returns a Selector that is true for particles from the hard vertex.
/// Particles without PU14 user info (e.g. ghost particles) will fail
/// the test.
inline 
fastjet::Selector SelectorIsHard()   {return SelectorVertexNumber(0);}
/// alternative name for SelectorHard
inline 
fastjet::Selector SelectorHard()   {return SelectorIsHard();}

/// returns a Selector that is true for particles from pileup vertices
/// and for particles without PU14 user info
inline 
fastjet::Selector SelectorIsPileup() {return !SelectorVertexNumber(0);}
/// alternative name for SelectorIsPileup
inline 
fastjet::Selector SelectorPileup() {return SelectorIsPileup();}

/// returns a Selector that is true for particles with a specific PDGId
fastjet::Selector SelectorPDGId(int i);

/// returns a Selector that is true for particles with a specific absolute PDGId
fastjet::Selector SelectorAbsPDGId(int i);


#endif
