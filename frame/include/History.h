#ifndef HISTORY_H
#define HISTORY_H

#include <sstream>
#include <iostream>
#include <string>

/// History class.
/// The class stores the history of processing steps.
/// New texts are added as strings. If setVerbosity(1) has been called, 
/// all the texts are displayed at command line when calling add().

class History {
 public:
  /// Default constructor.
  History();
  /// Argumented constructor.
  History(std::string text);
  /// Copy constructor.
  History(const History&);
  /// Copy operator.
  void operator= (const History&);
  /// Clear the History.
  void clear();
  /// Append a text to the History.
  void append(std::string text);
  /// Append a text to the History.
  /// After appending text to History, text gets cleared.
  void append(std::ostringstream& text);
  /// Set the verbosity of add().
  /// If set to 1, all texts added to History are printed to stdout.
  /// The default value is 0.
  static void setVerbosity(bool verb);
  /// Get the actual verbosity.
  bool getVerbosity() const;
  /// Get the whole History string.
  std::string getContent() const;
  
private:
  std::ostringstream history;
  static bool verbose;
};

#endif
