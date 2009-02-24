#include <ShapeLensConfig.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>
#include <limits>

using namespace std;

bool ShapeLensConfig::VERBOSITY = 0;
unsigned int ShapeLensConfig::NMAX_LOW = 0;
unsigned int ShapeLensConfig::NMAX_HIGH = 100;
data_t ShapeLensConfig::BETA_LOW = 0;
data_t ShapeLensConfig::BETA_HIGH = numeric_limits<data_t>::infinity();
data_t ShapeLensConfig::DELTA_BETA = 0.02;
bool ShapeLensConfig::ALLOW_FLATTENING = 0;
bool ShapeLensConfig::FILTER_SPURIOUS = 0;
data_t ShapeLensConfig::ADD_BORDER = 0.5;
unsigned int ShapeLensConfig::MIN_PIXELS = 20;
data_t ShapeLensConfig::MIN_THRESHOLD = 1.25;
data_t ShapeLensConfig::DETECT_THRESHOLD = 3.;
std::string ShapeLensConfig::NOISEMODEL = "GAUSSIAN";
bool ShapeLensConfig::BLENDING = 1;
data_t ShapeLensConfig::BLEND_MINCONT = 0.01;
unsigned int ShapeLensConfig::BLEND_NTHRESH = 16;
bool ShapeLensConfig::PIXEL_INTEGRATION = 0;

ShapeLensConfig::ShapeLensConfig() {
}

ShapeLensConfig::ShapeLensConfig(string filename) {
  // open config file
  ifstream configfile (filename.c_str());
  if (configfile.fail()) {
    cerr << "ShapeLensConfig: configuration file does not exists!" << endl;
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
	VERBOSITY = boost::lexical_cast<bool>(column[1].c_str());
      if (column[0] == "NMAX_LOW")
	NMAX_LOW = boost::lexical_cast<unsigned int>(column[1].c_str());
      if (column[0] == "NMAX_HIGH")
	NMAX_HIGH = boost::lexical_cast<unsigned int>(column[1].c_str());
      if (column[0] == "BETA_LOW")
	BETA_LOW = boost::lexical_cast<data_t>(column[1].c_str());
      if (column[0] == "BETA_HIGH")
	BETA_HIGH = boost::lexical_cast<data_t>(column[1].c_str());
      if (column[0] == "DELTA_BETA")
	DELTA_BETA = boost::lexical_cast<data_t>(column[1].c_str());
      if (column[0] == "ALLOW_FLATTENING")
	ALLOW_FLATTENING = boost::lexical_cast<bool>(column[1].c_str());
      if (column[0] == "FILTER_SPURIOUS")
        FILTER_SPURIOUS = boost::lexical_cast<bool>(column[1].c_str());
      if (column[0] == "ADD_BORDER")
        ADD_BORDER = boost::lexical_cast<data_t>(column[1].c_str());
      if (column[0] == "MIN_PIXELS")
	MIN_PIXELS = boost::lexical_cast<unsigned int>(column[1].c_str());
      if (column[0] == "MIN_THRESHOLD")
	MIN_THRESHOLD = boost::lexical_cast<data_t>(column[1].c_str());
      if (column[0] == "DETECT_THRESHOLD")
	DETECT_THRESHOLD = boost::lexical_cast<data_t>(column[1].c_str());
      if (column[0] == "NOISEMODEL")
	NOISEMODEL = column[1];
      if (column[0] == "PIXEL_INTEGRATION")
	PIXEL_INTEGRATION = boost::lexical_cast<bool>(column[1].c_str());
    }
  }
  // check if limits on beta and nmax make sense
  if (BETA_HIGH < BETA_LOW) {
    data_t tmp = BETA_LOW;
    BETA_LOW = BETA_HIGH;
    BETA_HIGH = tmp;
  }
  if (NMAX_HIGH < NMAX_LOW) {
    unsigned int tmp = NMAX_LOW;
    NMAX_LOW = NMAX_HIGH;
    NMAX_HIGH = tmp;
  }
}

