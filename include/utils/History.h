#ifndef HISTORY_H
#define HISTORY_H

#include <sstream>
#include <iostream>
#include <string>
#include <ShapeLensConfig.h>

namespace shapelens {

/// History class.
/// The class stores the history of processing steps.
/// New texts are added by using the operator<< (borrowed from std::ostringstream).\n
/// If setVerbosity(1) has been called, everything appended to History is also printed to stdout.

class History {
 public:
  /// Default constructor.
  History() {
    s.clear();
    silent = 0;
  }
  /// Copy constructor
  History(const History& h) {
    s.clear();
    s << h.s.str();
    silent = h.silent;
  }
  /// Copy operator.
  void operator=(const History& h) {
    s.clear();
    s << h.s.str();
    silent = h.silent;
  }
  /// Overloaded operator<<.
  /// With this operator, History behaves like a std::ostringstream.
  template <typename T> History& operator<<(T t) {
    if (ShapeLensConfig::VERBOSITY && !silent)
      std::cout << t;
    s << t; 
    return *this; 
  }

  History& operator<<(std::ostream& (*func)(std::ostream&)) {
    if (ShapeLensConfig::VERBOSITY && !silent)
      std::cout << func;
    s << func;
    return *this;
  } 

  History& operator<<(std::ios& (*func)(std::ios&)) {
    if (ShapeLensConfig::VERBOSITY && !silent)
      std::cout << func;
    s << func;
    return *this;
  }

  History& operator<<(std::ios_base& (*func)(std::ios_base&)) {
    if (ShapeLensConfig::VERBOSITY && !silent)
      std::cout << func;
    s << func;
    return *this;
  }
  /// Return the history string.
  std::string str() const {
    return s.str();
  }
  /// Clear the history.
  void clear() {
    s.str("");
  }
  /// Set the verbosity of operator<<.
  /// If set to 1, all texts added to History are printed to stdout.
  /// The default value is 0.
  static void setVerbosity(bool v) {
    ShapeLensConfig::VERBOSITY = v;
  }
  /// Get the actual verbosity.
  bool getVerbosity() const {
    return ShapeLensConfig::VERBOSITY;
  }
  /// Turn verbose mode off temporarily.
  /// This does not change <tt>verbosity</tt>, but turns
  /// the output off in cases where it is unwanted.
  void setSilent() {
    silent = 1;
  }
  /// Undo setSilent().
  void unsetSilent() {
    silent = 0;
  }
  /// Wheter history is empty.
  bool isEmpty() const {
    return !(bool) s.str().size();
  }
 private:
  std::ostringstream s;
  bool silent;
}; 
} // end namespace
#endif
