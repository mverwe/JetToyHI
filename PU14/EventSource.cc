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
                                    double &event_weight, double &prodX, double &prodY,
                                    int vertex_number) {
  PseudoJet particle;
  string line;
  double px, py, pz, m, E;
  int pdgid, vertex;

  unsigned original_size = particles.size();
  event_weight = 1;

  // read in particles
  while (getline(*_stream, line)) {
    // ignore blank lines and comment lines
    if (line.length() == 0 || line[0] == '#') continue;

    // if the line says "end" then assume we've found the end of the
    // event (try to make the check as efficient as possible).
    if (line[0] == 'e' && line.substr(0,3) == "end") break;

    bool bStop = false;
    // if the line says "weight", we multiply current weight by the number that follows
    if(line[0] == 'w' && line.substr(0, 6) == "weight")
    {
       string dummy;
       double temp = 1;
       FastIStringStream readline(line.c_str());
       readline >> dummy >> temp;
       if(temp > 0)
          event_weight = event_weight * temp;
       bStop = true;
    }

    //extraction production point info for hybrid event generator
    if(line.find("cross") < line.size()) {
      std::size_t posX = line.find("X");
      std::string strX = line.substr(posX+2,6);
      
      std::size_t posY = line.find("Y");
      std::string strY = line.substr(posY+2,6);
      
      prodX = atof(strX.c_str());
      prodY = atof(strY.c_str());

      bStop = true;
    }

    if(bStop) continue;

    // FastIStringStream is not a proper stream, but it's a lot faster
    // than standard istringstream.
    FastIStringStream readline(line.c_str());
    readline >> px >> py >> pz >> m >> pdgid >> vertex;
    assert(!readline.error());
    
    E = sqrt(px*px + py*py + pz*pz + m*m);
    particle = PseudoJet(px,py,pz,E);

    // now set the user info
    int barcode = particles.size();
    particle.set_user_info(new PU14(pdgid, barcode, vertex));

    // and add the particle to our final output
    particles.push_back(particle);
  }

  // if there were no new particles, then we assume the event has an error
  return (particles.size() != original_size);
}

