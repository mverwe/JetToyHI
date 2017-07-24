#ifndef __EVENTMIXER_HH__
#define __EVENTMIXER_HH__

#include "CmdLine.hh"
#include "EventSource.hh"
#include "GSLRandom.hh"

//----------------------------------------------------------------------
/// \class EventMixer
/// 
/// Class for that produces merged events from a mix of individual
/// hard and pileup interactions. It deduces the hard and pileup event
/// filenames from the command line, with options:
///
///   -hard HardFileName  -pileup PileupFileName  -npu NPU  [-chs]
///
/// Currently the number of pileup events is fixed (not Poisson
/// distributed).
///
/// Additional options:
///  -chs       when present, the charged pileup particles come scaled
///             by a factor 10^{-60} (this factor can be modified from
///             within code).
///  -massless  when present, particles come massless
///
///
/// Pileup multiplicities can be specified using one of the follosing
/// three options (by order of priority)
///  -npu  <N>   fixed npu=N
///  -mupu <N>   Poisson-distributed with average N
///  -upu  <N>   uniform between 1 and N
class EventMixer {
public:
  EventMixer(CmdLine * cmdline);
  ~EventMixer() {}

  /// causes the next event to be read in and mixed (hard + multiple pileup).
  /// Returns true if it successfully produced the event, false otherwise.
  bool next_event(); 

  /// returns a reference to vector of particles in the last event
  /// that was read in
  const std::vector<fastjet::PseudoJet> & particles() const {return _particles;}

  /// returns the number of pileup events generated in the last mixed event 
  int npu() const {return _npu;}

  /// Charged-hadron subtraction (CHS) can be "simulated" by scaling
  /// the charged particles from PU vertices by a factor
  /// _chs_rescaling_factor chosen << 1. This function returns the
  /// current value of that factor.
  double chs_rescaling_factor() const {return _chs_rescaling_factor;}

  /// set the chs_rescaling factor. If this is equal to one (default),
  /// then no rescaling is performed
  void set_chs_rescaling_factor(double r) {_chs_rescaling_factor = r;}

  /// return true if events are generated using massless particles
  bool massless() const { return _massless;}

  /// returns a description of what the EventMixer does (useful for bookkeeping)
  std::string description() const;

private:
  CmdLine * _cmdline;
  std::string _hard_name, _pileup_name;
  fastjet::SharedPtr<EventSource> _hard, _pileup;
  int _npu, _mupu, _upu;
  GSLRandom _rng;
  double _chs_rescaling_factor;
  bool _massless;

  std::vector<fastjet::PseudoJet> _particles;
};

#endif  // __EVENTMIXER_HH__

