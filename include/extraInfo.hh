#ifndef extraInfo_h
#define extraInfo_h

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

using namespace std;
using namespace fastjet;

class extraInfo : public PseudoJet::UserInfoBase{
public:
  // default ctor
  //  - pdg_id        the PDG id of the particle
  //  - vertex_number the id of the vertex it originates from
  extraInfo(const int & pdg_id_in, const int & vertex_number_in) :
    _pdg_id(pdg_id_in), _vertex_number(vertex_number_in){}

  /// access to the PDG id
  int pdg_id() const { return _pdg_id;}

  /// access to the vertex number
  int vertex_number() const { return _vertex_number;}

protected:
  int _pdg_id;         // the associated pdg id
  int _vertex_number;  // the associated vertex number
};

#endif
