#include <History.h>

History::History() {
  clear();
}

History::History(std::string text) {
  clear();
  history << text;
}

History::History(const History& H) {
  clear();
  history << H.getContent();
  verbose = H.getVerbosity();
}

void History::operator= (const History& H) {
  clear();
  history << H.getContent();
  verbose = H.getVerbosity();
}

bool History::verbose=0;

void History::clear() {
  history.str("");
}

void History::append(std::string text) {
  history << text;
  if (verbose)
    std::cout << text;
}

void History::append(std::ostringstream& text) {
  append(text.str());
  text.str("");
}

void History::setVerbosity(bool verb) {
  verbose = verb;
}

bool History::getVerbosity() const {
  return verbose;
}

std::string History::getContent() const {
  return history.str();
}

