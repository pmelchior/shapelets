#include <ShapeLensConfig.h>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <iostream>

using namespace std;

unsigned int ShapeLensConfig::NMAX_LOW = 0;
unsigned int ShapeLensConfig::NMAX_HIGH = 100;
double ShapeLensConfig::BETA_LOW = 0;
double ShapeLensConfig::BETA_HIGH = INFINITY;
bool ShapeLensConfig::REGULARIZE = 0;
double ShapeLensConfig::REG_LIMIT = 1e-5;
std::string ShapeLensConfig::UNREG_SIFFILE = "";
bool ShapeLensConfig::ALLOW_FLATTENING = 0;
bool ShapeLensConfig::FILTER_SPURIOUS = 0;
double ShapeLensConfig::ADD_BORDER = 0.5;

ShapeLensConfig::ShapeLensConfig() {
}

ShapeLensConfig::ShapeLensConfig(string filename) {
  // open config file
  ifstream configfile (filename.c_str());
  if (configfile.fail()) {
    cout << "ShapeLensConfig: configuration file does not exists!" << endl;
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
    // exclude empty and comment lines
    if (column.size() >= 2 && column[0].compare("#") != 0) {
      if (column[0].compare("NMAX_LOW") == 0)
	NMAX_LOW = (unsigned int) atoi (column[1].c_str());
      if (column[0].compare("NMAX_HIGH") == 0)
	NMAX_HIGH = (unsigned int) atoi (column[1].c_str());
      if (column[0].compare("BETA_LOW") == 0)
	BETA_LOW = (double) atof (column[1].c_str());
      if (column[0].compare("BETA_HIGH") == 0)
	BETA_HIGH = (double) atof (column[1].c_str());
      if (column[0].compare("REGULARIZE") == 0)
	REGULARIZE = (bool) atoi (column[1].c_str());
      if (column[0].compare("REG_LIMIT") == 0)
	REG_LIMIT = (double) atof (column[1].c_str());
      if (column[0].compare("UNREG_SIFFILE") == 0)
	UNREG_SIFFILE = column[1];
      if (column[0].compare("ALLOW_FLATTENING") == 0)
	ALLOW_FLATTENING = (bool) atoi (column[1].c_str());
      if (column[0].compare("FILTER_SPURIOUS") == 0)
        FILTER_SPURIOUS = (bool) atoi (column[1].c_str());
      if (column[0].compare("ADD_BORDER") == 0)
        ADD_BORDER = (double) atof (column[1].c_str());
    }
  }
}

