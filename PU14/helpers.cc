#include "helpers.hh"

using namespace fastjet;
using namespace std;

double Width::result(const fastjet::PseudoJet &jet) const{
  // check the jet is appropriate for computation
  if (!jet.has_constituents())
    throw Error("Angularities can only be applied on jets for which the constituents are known.");
  
  vector<PseudoJet> constits = jet.constituents();
  double num=0.0, den=0.0;
  for (vector<PseudoJet>::iterator ci = constits.begin(); ci!=constits.end(); ci
++){
    double pt = ci->pt();
    num += pt * ci->delta_R(jet);
    den += pt;
  }
  return num/den;

}
