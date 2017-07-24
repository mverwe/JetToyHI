#ifndef __GSLRANDOM_HH__
#define __GSLRANDOM_HH__

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#ifdef USE_CXX11
#include <memory>
#define SHARED_PTR std::shared_ptr
#else
#include <tr1/memory>
#define SHARED_PTR std::tr1::shared_ptr
#endif

#include <string>

/// A (partial) C++ interface to the GSL random number generators;
///
/// - to link remember to use -lgsl -lgslcblas
///
/// - names of generators are to be found at
///     http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html
///   & include gsl_rng_mt19937 (default) 
///             gsl_rng_ranlxs0 / s1 / s2 (second generation ranlux)
///             gsl_rng_ranlxd1 / d2      (48 bits, double precision output)
///             gsl_rng_ranlux / ranlux289 (orig. ranlux; 389: 24 decorr bits)
///             gsl_rng_cmrg (combined multiple recursive gen [l'Ecuyer])
///      etc.... (see docs for more info)
///
/// - some enquiry functions haven't been implemented (e.g. for state, range).
///
class GSLRandom {
public:
  inline GSLRandom() {
    reset(gsl_rng_alloc(gsl_rng_mt19937));
  }
  inline GSLRandom(unsigned long int s) {
    reset(gsl_rng_alloc(gsl_rng_mt19937));
    set(s);
  }
  inline GSLRandom(const gsl_rng_type * T) {reset(gsl_rng_alloc(T));}
  inline GSLRandom(const gsl_rng_type * T, unsigned long int s) {
    reset(gsl_rng_alloc(T));
    set(s);
  }

  inline void reset(gsl_rng * r) {_r.reset(r, gsl_rng_free);}

  /// set the seed
  inline void set(unsigned long int s) {gsl_rng_set(_r.get(), s);}

  /// returns in range [0,1) (includes 0, excludes 1)
  inline double uniform() const {return gsl_rng_uniform(_r.get());}

  /// returns in range (0,1) (excludes 0, excludes 1)
  inline double uniform_pos() const {return gsl_rng_uniform_pos(_r.get());}

  /// returns in range [xmin,xmax) (includes xmin, excludes xmax -- modulo rounding)
  inline double uniform(double xmin, double xmax) const {return xmin + (xmax-xmin)*gsl_rng_uniform(_r.get());}

  /// returns +-1 with equal probability
  inline double sign() const {return uniform() < 0.5 ? -1 : 1;}

  /// returns a Gaussian distributed random number with standard
  /// deviation sigma and mean zero
  inline double gaussian(double sigma) const {return gsl_ran_gaussian(_r.get(),sigma);}
  /// Gaussian with sigma=1
  inline double gaussian() const {return gsl_ran_ugaussian(_r.get());}

  /// return random x according to p(x) dx = {1 \over \mu} \exp(-x/\mu) dx
  inline double exponential(double mu) const {return gsl_ran_exponential(_r.get(),mu);}

  // returns random k according to p(k) = {\mu^k \over k!} \exp(-\mu)
  inline unsigned int poisson(double mu) const {return gsl_ran_poisson(_r.get(),mu);}

  // returns a number generated according to the gamma distribution
  inline double gamma(double a, double b) const { return gsl_ran_gamma(_r.get(),a,b);}

  /// returns integer in range [0,n-1]
  inline unsigned long int uniform_int(unsigned long int n) const {
    return n > 1 ? gsl_rng_uniform_int(_r.get(), n) : 0;}

  inline std::string name() const {return std::string(gsl_rng_name (_r.get()));}

  inline gsl_rng * gsl_generator() {return _r.get();}
  
  /// no specific destructor
  //inline ~GSLRandom() {}
protected:
  SHARED_PTR<gsl_rng> _r;
};

#endif // __GSLRANDOM_HH__
