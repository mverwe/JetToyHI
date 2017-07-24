// -*- C++ -*-
#ifndef __FASTISTRINGSTREAM_HH__
#define __FASTISTRINGSTREAM_HH__

#include <string>
#include <cstdlib>

/// class to help with conversion of the contents of a character array
/// into numeric data. It provides a significantly faster, (but much
/// more limited) alternative to the use of an istringstream.
///
/// Its implementation is inspired by the discussion at
/// http://stackoverflow.com/questions/5678932/fastest-way-to-read-numerical-values-from-text-file-in-c-double-in-this-case
///
/// Note that this class does _not_ currently derive from istream.
/// 
/// Copyright (C) 2012-2013 Gavin P. Salam
/// Released under the terms of the GNU General Public License v3
class FastIStringStream {
public:
  /// constructor 
  FastIStringStream(const char * line) : _next(const_cast<char *>(line)), 
                                _new_next(_next), 
                                _error(false) {}

  // while this is templated, it will currently only work with 
  // double, float and int
  template<class T> FastIStringStream & operator>>(T & value) {
    _get(value);
    // the following line exploits the fact that if strtod and strtol,
    // below, fail to read a number then they leave _new_next equal to
    // _next (as long as _new_next was not zero to start with).
    if (_new_next == _next) _error = true;
    _next = _new_next;
    return *this;
  }

  /// returns true if an error was encountered during any of the
  /// operator>> calls
  bool error() {return _error;}

  /// @brief Quick status checking
  /// @return true on absence of any errors
  ///
  /// this operator allows the class to be interpreted as a boolean,
  /// e.g. within an if statement to see if a read was successful
  operator bool() const {return !_error;}

private:
  void _get(double & x) {
    x = std::strtod(_next, &_new_next);
  }
  void _get(float & x) {
    x = std::strtof(_next, &_new_next);
  }
  void _get(int & i) {
    i = std::strtol(_next, &_new_next, 10); // force base 10
  }
  void _get(std::string & s) {
    s.clear();
    bool started = false;
    while (*_new_next != 0) {
      if (*_new_next == ' ') {
        if (started) break;
        // otherwise simply skip leading space
      } else {
        started = true;
        s.append(1, *_new_next);
      }
      ++_new_next;
    }
  }

  /// a pointer to the next character
  char * _next;
  /// after reading a number, the pointer to the next new character
  char * _new_next;
  bool _error;
};


#endif // __FASTISTRINGSTREAM_HH__
