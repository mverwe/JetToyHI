#ifndef __Angularity_HH__
#define __Angularity_HH__

//------------------------------------------------------------------------
/// angularities
///
/// This is defined as in https://arxiv.org/pdf/1408.3122.pdf
/// \f[
///   \lambda^{\kappa}_{\beta} = \sum_i (p_{ti}/p_{T,jet})^{\kappa} (\Delta R_{i,jet}/R_{0})^{\beta}
/// \f]
class Angularity {
public:
  /// default ctor
  Angularity(double beta=1.0, double kappa = 1., double R0 = 0.4) :
    _beta(beta),
    _kappa(kappa),
    _R0(R0)
  {}

  /// compute the function
  virtual double result(const fastjet::PseudoJet &jet) const {
    // check the jet is appropriate for computation
    if (!jet.has_constituents()) {
      Printf("Angularities can only be applied on jets for which the constituents are known.");
      return -999.;
    }
    vector<fastjet::PseudoJet> constits = jet.constituents();
    double ang = 0.;
    for(fastjet::PseudoJet p : constits) {
      ang += std::pow(p.pt()/jet.pt(),_kappa)*std::pow(jet.delta_R(p)/_R0,_beta); //use scalar pT for jet?
    }
    return ang;
  }
  
protected:
  double _beta;
  double _kappa;
  double _R0;
};

#endif
