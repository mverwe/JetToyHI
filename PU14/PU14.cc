#include "PU14.hh"
#include <sstream>
#include <memory>

using namespace std;
using namespace fastjet;

std::ostream & operator<<(std::ostream & o, const fastjet::PseudoJet & p) {
  o << "pt = " << p.pt() << ", "
    << "rap = " << p.rap() << ", "
    << "phi = " << p.phi() << ", "
    << "m = " << p.m() ;
  if (p.has_user_info<PU14>()) {
    o << ", pdg_id = " << p.user_info<PU14>().pdg_id();
    o << ", vertex = " << p.user_info<PU14>().vertex();
  }
  return o;
}


//----------------------------------------------------------------------
/// the worker class for charged particle selection
class SelectorWorkerIsCharged : public SelectorWorker {
public:

  virtual bool pass(const PseudoJet & particle) const {
    // we check that the user_info_ptr is non-zero so as to make
    // sure that explicit ghosts don't cause the selector to fail
    return (particle.user_info_ptr() != 0 && 
	    particle.user_info<PU14>().three_charge() != 0);
  }
  
  virtual string description() const {return "is_charged";}
};

//----------------------------------------------------------------------
Selector SelectorIsCharged() {
  return new SelectorWorkerIsCharged();
}


//----------------------------------------------------------------------
/// the worker class for particle selection according to the internal
/// vertex number
class SelectorWorkerVertexNumber : public SelectorWorker {
public:
  SelectorWorkerVertexNumber(int vertex_number) : _vertex_number(vertex_number) {}

  virtual bool pass(const PseudoJet & particle) const {
    // we check that the user_info_ptr is non-zero so as to make
    // sure that explicit ghosts don't cause the selector to fail
    return (particle.user_info_ptr() != 0 && 
	    particle.user_info<PU14>().vertex_number() == _vertex_number);
  }
  
  virtual string description() const {
    ostringstream ostr;
    ostr << "vertex number == " << _vertex_number;
    return ostr.str();
  }

private:
  int _vertex_number;
};

Selector SelectorVertexNumber(int vertex_number) {
  return new SelectorWorkerVertexNumber(vertex_number);
}



// //----------------------------------------------------------------------
// /// the worker class for particles that have net B-flavour
// class SelectorWorkerHasBFlavour : public SelectorWorker {
// public:
// 
//   virtual bool pass(const PseudoJet & particle) const {
//     // we check that the user_info_ptr is non-zero so as to make
//     // sure that explicit ghosts don't cause the selector to fail
//     return (particle.user_info_ptr() != 0 && 
// 	    particle.user_info<ParticleInfo>().flav[5] != 0);
//   }
//   
//   virtual string description() const {return "has_b_flavour";}
// };
// 
// 
// //----------------------------------------------------------------------
// Selector SelectorHasBFlavour() {
//   return new SelectorWorkerHasBFlavour();
// }


//----------------------------------------------------------------------
class SelectorWorkerPDGId : public SelectorWorker {
public:
  SelectorWorkerPDGId(int i) : _id(i) {}
  virtual bool pass(const PseudoJet & particle) const {
    return (particle.user_info_ptr() != 0 && 
	    particle.user_info<PU14>().pdg_id() == _id);
  }
  virtual string description() const {
    ostringstream ostr;
    ostr << "pdg_id==" << _id;
    return ostr.str();
  }
private:
  int _id;
};

Selector SelectorPDGId(int i) {return new SelectorWorkerPDGId(i);}

//----------------------------------------------------------------------
class SelectorWorkerAbsPDGId : public SelectorWorker {
public:
  SelectorWorkerAbsPDGId(int i) : _id(i) {}
  virtual bool pass(const PseudoJet & particle) const {
    return (particle.user_info_ptr() != 0 && 
	    abs(particle.user_info<PU14>().pdg_id()) == _id);
  }
  virtual string description() const {
    ostringstream ostr;
    ostr << "|pdg_id|==" << _id;
    return ostr.str();
  }
private:
  int _id;
};

Selector SelectorAbsPDGId(int i) {return new SelectorWorkerAbsPDGId(i);}


// //----------------------------------------------------------------------
// class SelectorWorkerComesFrom : public SelectorWorker {
// public:
//   SelectorWorkerComesFrom(HepMC::GenEvent * gen_event, const Selector & from) {
//     _from = from;
//     // first get the list of particles that satisfy the from condition
//     for (HepMC::GenEvent::particle_iterator pit = gen_event->particles_begin();
//          pit != gen_event->particles_end();
//          pit++) {
//       PseudoJet p(**pit);
//       if (from.pass(p)) _from_particles.push_back(*pit);
//     }
//     // and now create an object that essentially contains a set of all
//     // those particles' descendants to decide subsequently, whether
//     // any subsequently supplied particle forms part of that set
//     _is_x_not_from.reset(new IsXNotFromY<IsAlwaysTrue>(IsAlwaysTrue(), _from_particles));
//   }
//   virtual bool pass(const PseudoJet & particle) const {
//     if (! particle.has_user_info<PI>()) return false;
//     const HepMC::GenParticle * gen_p = particle.user_info<PI>().hepmc_particle;
//     if (! gen_p) return false;
//     return ! (*_is_x_not_from)(gen_p);
//   }
//   virtual string description() const {
//     ostringstream ostr;
//     ostr << "particles descending from the " 
//          << _from_particles.size() 
//          << " that satisfied (" 
//          << _from.description() << ")";
//     return ostr.str();
//   }
// private:
//   vector<HepMC::GenParticle *> _from_particles;
//   std::auto_ptr<IsXNotFromY<IsAlwaysTrue> > _is_x_not_from;
//   Selector _from;
// };
// 
// Selector SelectorComesFrom(HepMC::GenEvent * gen_event, 
//                            const fastjet::Selector & from) {
//   return new SelectorWorkerComesFrom(gen_event, from);
// }
