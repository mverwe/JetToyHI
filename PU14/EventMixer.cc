#include "EventMixer.hh"
#include "PU14.hh"
#include "helpers.hh"
#include <cstdlib>

using namespace std;

//----------------------------------------------------------------------
EventMixer::EventMixer(CmdLine * cmdline) : _cmdline(cmdline) {

  _hard_name = _cmdline->value<string>("-hard");
  _pileup_name = _cmdline->value<string>("-pileup", "");

  // 3 ways of setting the multiplicity of pileup events
  //
  //  -npu <npu>  : fixed <npu> number of PU vertices
  //  -upu <upu>  : uniformly-distributed number of PU
  //                vertices between 1 and upu
  //  -mupu <mu>  : Poisson-distributed (with average mu)
  //                number of pileup vertices

  // uniform
  _upu = _cmdline->value("-upu", -1);

  // poisson-distriuted
  // Note that _mupu could in principle be a double
  if (_upu < 0){
    _mupu = _cmdline->value("-mupu", -1);
  } else {
    _mupu = -1;
  }

  // fixed (the default)
  if ((_upu < 0) && (_mupu < 0)){
    // only allow the -npu option if _upu is not set
    _npu = _cmdline->value("-npu", 20);
  } else {
    _npu = -1; // will be reset by event mixer
  }
  
  
  _massless = _cmdline->present("-massless");

  if (_cmdline->present("-chs")) {
    set_chs_rescaling_factor(1e-60);
  } else {
    set_chs_rescaling_factor(1.0); // this effectively turns off CHS
  }

  _hard  .reset(new EventSource(_hard_name  ));

  if (_pileup_name.empty()){
    cerr << "WARNING: no pileup included" << endl;
    _pileup.reset();
    _mupu=_upu=0;
    _npu=0;
  } else {
    _pileup.reset(new EventSource(_pileup_name));
  }

  // initialise a possible random number generator
  _rng.set(_cmdline->value<unsigned long int>("-seed", 1));
}

//----------------------------------------------------------------------
bool EventMixer::next_event() {
  _particles.resize(0);
  
  // first get the hard event
  if (! _hard->append_next_event(_particles,0)) return false;

  unsigned hard_size = _particles.size();

  // add pileup if available
  if (_pileup.get()){
    // allow for random
    // Uniform case:
    if (_upu > 0) {
      //_npu = 1 + int(_upu*(1.0*rand()/RAND_MAX));
      _npu = 1+_rng.uniform_int(_upu);
      if (_npu > _upu) _npu = _upu; // for the case of equality with RAND_MAX
    }
    // Poisson-distributed
    if (_mupu>0){
      _npu = _rng.poisson(_mupu); 
    }      
    for (int i = 1; i <= _npu; i++) {
      if (! _pileup->append_next_event(_particles,i)) return false;
    }
  }

  // make particles massless if requested
  if (_massless){
    _particles = MasslessTransformer()(_particles);
  }

  // apply CHS rescaling factor if requested
  if (chs_rescaling_factor() != 1.0) {
    for (unsigned i = hard_size; i < _particles.size(); i++) {
      if (_particles[i].user_info<PU14>().charge() != 0) _particles[i] *= _chs_rescaling_factor;
    }
  }

  return true;
}


//----------------------------------------------------------------------
string EventMixer::description() const {
  ostringstream ostr;
  ostr << "Event mixer using hard events from " << _hard_name;

  if (_upu < 0) {
    ostr << " and " << _npu << " pileup events from " << _pileup_name;
  } else {
    ostr << " and 1..." << _upu << " pileup events from " << _pileup_name
         << " (exact number chosen according to uniform random distribution)";
  }
  if (chs_rescaling_factor() != 1.0) {
    ostr << " with CHS (rescaling charged PU by a factor "
         << chs_rescaling_factor() << ")";
  }
  if (_massless){
    ostr << " and massless particles";
  }
  return ostr.str();
}
