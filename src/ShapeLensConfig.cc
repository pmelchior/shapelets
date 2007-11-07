#include <ShapeLensConfig.h>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <iostream>
#include <limits>

using namespace std;

bool ShapeLensConfig::VERBOSITY = 0;
unsigned int ShapeLensConfig::NMAX_LOW = 0;
unsigned int ShapeLensConfig::NMAX_HIGH = 100;
data_t ShapeLensConfig::BETA_LOW = 0;
data_t ShapeLensConfig::BETA_HIGH = numeric_limits<data_t>::infinity();
bool ShapeLensConfig::REGULARIZE = 0;
data_t ShapeLensConfig::REG_LIMIT = 1e-5;
bool ShapeLensConfig::SAVE_UNREG = 1;
bool ShapeLensConfig::ALLOW_FLATTENING = 0;
bool ShapeLensConfig::FILTER_SPURIOUS = 0;
data_t ShapeLensConfig::ADD_BORDER = 0.5;
unsigned int ShapeLensConfig::MIN_PIXELS = 20;
data_t ShapeLensConfig::MIN_THRESHOLD = 1.25;
data_t ShapeLensConfig::DETECT_THRESHOLD = 3.;
std::string ShapeLensConfig::NOISEMODEL = "GAUSSIAN";

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
    if (column.size() >= 2 && column[0] != "#") {
      if (column[0] == "VERBOSITY")
	VERBOSITY = (bool) atoi (column[1].c_str());
      if (column[0] == "NMAX_LOW")
	NMAX_LOW = (unsigned int) atoi (column[1].c_str());
      if (column[0] == "NMAX_HIGH")
	NMAX_HIGH = (unsigned int) atoi (column[1].c_str());
      if (column[0] == "BETA_LOW")
	BETA_LOW = (data_t) atof (column[1].c_str());
      if (column[0] == "BETA_HIGH")
	BETA_HIGH = (data_t) atof (column[1].c_str());
      if (column[0] == "REGULARIZE")
	REGULARIZE = (bool) atoi (column[1].c_str());
      if (column[0] == "REG_LIMIT")
	REG_LIMIT = (data_t) atof (column[1].c_str());
      if (column[0] == "SAVE_UNREG")
	SAVE_UNREG = (bool) atoi(column[1].c_str());
      if (column[0] == "ALLOW_FLATTENING")
	ALLOW_FLATTENING = (bool) atoi (column[1].c_str());
      if (column[0] == "FILTER_SPURIOUS")
        FILTER_SPURIOUS = (bool) atoi (column[1].c_str());
      if (column[0] == "ADD_BORDER")
        ADD_BORDER = (data_t) atof (column[1].c_str());
      if (column[0] == "MIN_PIXELS")
	MIN_PIXELS = (unsigned int) atoi (column[1].c_str());
      if (column[0] == "MIN_THRESHOLD")
	MIN_THRESHOLD = (data_t) atof (column[1].c_str());
      if (column[0] == "DETECT_THRESHOLD")
	DETECT_THRESHOLD = (data_t) atof (column[1].c_str());
      if (column[0] == "NOISEMODEL")
	NOISEMODEL = column[1];
    }
  }
}

