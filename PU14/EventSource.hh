#ifndef __EVENTSOURCE_HH__
#define __EVENTSOURCE_HH__

#include <string>
#include <vector>
#include <istream>
#include <memory>
#include "fastjet/PseudoJet.hh"
#include "fastjet/SharedPtr.hh"

//----------------------------------------------------------------------
/// \class EventSource
///
/// Class for reading events from a file (or stdin) in some simple format
class EventSource {
public:
  
  EventSource(const std::string & filename) {
    open_stream(filename);
  }

  /// set up an event stream from the corresponding file (in the PU14 format) 
  void open_stream(const std::string & filename);

  /// appends the particles from the next event that is read onto the 
  /// particles vector.
  bool append_next_event(std::vector<fastjet::PseudoJet> & particles,
                         double &event_weight, double &prodX, double &prodY,
                         int vertex_number = 0);

private:
  std::istream * _stream;
  fastjet::SharedPtr<std::istream> _stream_auto;

};
 
#endif // __EVENTSOURCE_HH__
