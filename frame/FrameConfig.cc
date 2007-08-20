#include <FrameConfig.h>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <iostream>

using namespace std;

bool FrameConfig::FILTER_SPURIOUS = 0;
double FrameConfig::ADD_BORDER = 0.5;

FrameConfig::FrameConfig() {
}

FrameConfig::FrameConfig(string filename) {
  // open config file
  ifstream configfile (filename.c_str());
  if (configfile.fail()) {
    cout << "FrameConfig: configuration file does not exists!" << endl;
    terminate();
  }
  // read in config file
  string line;
  while(getline(configfile, line)) {
    typedef boost::tokenizer<boost::char_separator<char> > Tok;
    // split entries at tabs
    boost::char_separator<char> sep("\t");
    Tok tok(line, sep);
    // first of all we copy the token into string vector
    // though this is not too sophisticated it does not hurt and is more 
    // convenient
    std::vector<std::string> column;
    for(Tok::iterator tok_iter = tok.begin(); tok_iter != tok.end(); ++tok_iter)
      column.push_back(*tok_iter);
    // exclude comment and empty lines
    if (column.size() >= 2 && column[0].compare("#") != 0) {
      if (column[0].compare("FILTER_SPURIOUS") == 0)
	FILTER_SPURIOUS = (bool) atoi (column[1].c_str());
      if (column[0].compare("ADD_BORDER") == 0)
	ADD_BORDER = (double) atof (column[1].c_str());
    }
  }
}
