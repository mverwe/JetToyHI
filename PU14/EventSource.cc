#include "EventSource.hh"
#include "FastIStringStream.hh"
#include "PU14.hh"
#include "zfstream.h"
#include <cassert>
#include <fstream>

using namespace std;
using namespace fastjet;


//----------------------------------------------------------------------
void EventSource::open_stream(const std::string & filename) {
  if (filename == "-") {
    _stream  = & cin;
  } else if (filename.length() > 3 && 
             filename.find(std::string(".gz")) +3 == filename.length()) {
    _stream = new gzifstream(filename.c_str());
    _stream_auto.reset(_stream);
  } else {
    _stream = new ifstream(filename.c_str());
    _stream_auto.reset(_stream);
  }
  if (! _stream->good()) {
    cerr << "ERROR: could not access file " << filename << endl;
    exit(-1);
  }
}


//----------------------------------------------------------------------
bool EventSource::append_next_event(std::vector<fastjet::PseudoJet> & particles,
                                    int vertex_number) {
  PseudoJet particle;
  string line;
  double px, py, pz, m, E;
  int pdgid;

  unsigned original_size = particles.size();

  // read in particles
  while (getline(*_stream, line)) {
    // ignore blank lines and comment lines
    if (line.length() == 0 || line[0] == '#') continue;

    // if the line says "end" then assume we've found the end of the
    // event (try to make the check as efficient as possible).
    if (line[0] == 'e' && line.substr(0,3) == "end") break;

    // FastIStringStream is not a proper stream, but it's a lot faster
    // than standard istringstream.
    FastIStringStream readline(line.c_str());
    readline >> px >> py >> pz >> m >> pdgid;
    assert(!readline.error());
    
    E = sqrt(px*px + py*py + pz*pz + m*m);
    particle = PseudoJet(px,py,pz,E);

    // now set the user info
    int barcode = particles.size();
    particle.set_user_info(new PU14(pdgid, barcode, vertex_number));

    // and add the particle to our final output
    particles.push_back(particle);
  }

  // if there were no new particles, then we assume the event has an error
  return (particles.size() != original_size);
}

