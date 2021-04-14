#ifndef __JetCharge_HH__
#define __JetCharge_HH__

//------------------------------------------------------------------------
/// jet charge
///
/// This is defined as in 

class JetCharge {
public:
  /// default ctor
  JetCharge(double kappa = 0.5, double ptmin = -1.) :
    _kappa(kappa),
    _ptmin(ptmin)
   {}

  /// compute the function
  virtual double result(const fastjet::PseudoJet &jet) const {
    // check the jet is appropriate for computation
    if (!jet.has_constituents()) {
      Printf("Angularities can only be applied on jets for which the constituents are known.");
      return -999.;
    }
    vector<fastjet::PseudoJet> constits = jet.constituents();
    double sumcharge = 0.;
    double sumpt = 0.;
    for(fastjet::PseudoJet p : constits) {
      if(p.perp()<_ptmin) continue;
      const int & ch = p.user_info<PU14>().charge(); //three_charge()
      //std::cout << "charge: " << ch << " three_charge: " << p.user_info<PU14>().three_charge() << "  pdg: " << p.user_info<PU14>().pdg_id() << std::endl;
      sumcharge += ch*std::pow(p.pt(),_kappa);
      sumpt += std::pow(p.pt(),_kappa);
    }
    double charge = sumcharge/sumpt;
    return charge;
  }
  
protected:
  double _kappa;
  double _ptmin;
};

#endif
