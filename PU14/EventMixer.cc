#include "EventMixer.hh"
#include "PU14.hh"
#include "helpers.hh"
#include <cstdlib>

using namespace std;

//----------------------------------------------------------------------
EventMixer::EventMixer(CmdLine * cmdline) : _cmdline(cmdline) {

  _hard_name = _cmdline->value<string>("-hard");
  _pileup_name = _cmdline->value<string>("-pileup", "");

  // setting the multiplicity of pileup events (background HI)
  //
  //  -npu <npu>  : fixed <npu> number of PU vertices - default to 1

  // fixed (the default)
  _npu = _cmdline->value("-npu", 1);
  if(_npu > 1)
  {
     cerr << "WARNING: number of background event requested = " << _npu << endl;
     cerr << "   make sure you actually want that!" << endl;
  }

  _massless = _cmdline->present("-massless");

  if (_cmdline->present("-chs")) {
    set_chs_rescaling_factor(1e-60);
  } else {
    set_chs_rescaling_factor(1.0); // this effectively turns off CHS
  }

  _hard  .reset(new EventSource(_hard_name  ));

  if (_pileup_name.empty()){
    cerr << "INFO: no background requested" << endl;
    _pileup.reset();
    _npu=0;
  } else {
    _pileup.reset(new EventSource(_pileup_name));
  }
}

//----------------------------------------------------------------------
bool EventMixer::next_event() {
  _particles.resize(0);
  _hard_event_weight = 1;
  _pu_event_weight = 1;
  
  // first get the hard event
  if (! _hard->append_next_event(_particles,_hard_event_weight,_posX,_posY,0)) return false;

  unsigned hard_size = _particles.size();

  // add pileup if available
  if (_pileup.get()){
    for (int i = 1; i <= _npu; i++) {
      double dumX, dumY;
      if (! _pileup->append_next_event(_particles,_pu_event_weight,dumX,dumY,i)) return false;
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

  if (_npu > 0) {
    ostr << " and " << _npu << " pileup events from " << _pileup_name;
  
    if (chs_rescaling_factor() != 1.0) {
    ostr << " with CHS (rescaling charged PU by a factor "
         << chs_rescaling_factor() << ")";
    }
  }
  if (_massless){
    ostr << " and massless particles";
  }
  return ostr.str();
}
