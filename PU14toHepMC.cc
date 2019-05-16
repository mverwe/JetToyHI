#include <exception>
#include <sstream>
#include <fstream>
#include <vector>
#include "HepMC/IO_BaseClass.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/Units.h"

using namespace std;
using namespace HepMC;

int exitusage(int status) {
  cerr << "Usage: PU14toHepMC {outputfile} [{idA}={eA} {idB}={eB}] "
    "{inputfile} [{inputfile} ...]\n";
  return status;
}

int exiterror(string mess) {
  cerr << mess << endl;
  return 2;
}

bool convertEvent(istream & is, GenVertex & vx) {
  string line;
  int n = 2;
  while ( getline(is, line) ) {
    if ( line == "end" ) return true;
    istringstream iss(line);
    double px = 0.0, py = 0.0, pz = 0.0, m = 0.0;
    int id = 0, status = 0;
    if ( !( iss >> px >> py >> pz >> m >> id >> status ) )
      return false;
    GenParticle * p =
      new GenParticle(FourVector(px, py, pz, sqrt(px*px+py*py*pz*pz+m*m)), id, 1);
    p->suggest_barcode(status*1000000 + ++n);
    vx.add_particle_out(p);
  }

}
 
int main(int argc, char ** argv) {
 

  if ( argc < 3 ) return exitusage(1);

  int ida = 2212;
  int idb = 2212;
  double ea = 6500.0;
  double eb = 6500.0;

  string output = argv[1];
  int iarg = 2;
  string beam = argv[iarg];
  if ( beam.find("=") != string::npos ) {
    ida = stoi(beam.substr(0, beam.find("=")));
    ea = stod(beam.substr(beam.find("=") + 1));
    ++iarg;
  }
  if ( argc < iarg + 1 ) exitusage(1);
  beam = argv[iarg];
  if ( beam.find("=") != string::npos ) {
    idb = stoi(beam.substr(0, beam.find("=")));
    eb = stod(beam.substr(beam.find("=") + 1));
    ++iarg;
  }

  vector<string> inputs;
  for ( int i = iarg; i < argc; ++i ) inputs.push_back(argv[i]);

  HepMC::IO_GenEvent hepmcio(output, std::ios::out);
  int neve = 0;

  for ( string file : inputs ) {
    cout << "Reading " << file;
    ifstream is(file.c_str());
    string line;
    while ( getline(is, line) ) {
      if ( line.substr(0, 7) == "# event" ) {
        GenEvent * e = new GenEvent();
        GenVertex * vx = new GenVertex();
        GenParticle * ba =
          new GenParticle(FourVector(0.0, 0.0, ea, ea), ida, 4);
        GenParticle * bb =
          new GenParticle(FourVector(0.0, 0.0, -eb, eb), idb, 4);
        ba->suggest_barcode(1);
        bb->suggest_barcode(2);
        vx->add_particle_in(ba);
        vx->add_particle_in(bb);
        e->set_event_number(neve++);
        if ( !convertEvent(is, *vx) ) return exiterror("Failed to convert event");
        e->add_vertex(vx);
        e->set_beam_particles(ba, bb);
        e->set_signal_process_vertex(vx);
        if ( !(hepmcio << e) ) return exiterror("Failed to write event");
        delete e;
      }
    }
  }

  cout << "Wrote " << neve << " events to " << output << endl;

  return 0;
}

