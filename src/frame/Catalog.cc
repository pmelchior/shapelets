#include <frame/Catalog.h>
#include <fstream>
#include <iostream>
#include <boost/tokenizer.hpp>

using namespace std;

Catalog::Catalog() : vector<CatObject>() {
}

Catalog::Catalog(string catfile) : vector<CatObject>() {
  read(catfile);
}

void Catalog::read(string catfile) {
  // reset bitset that indicates which data are given in catalog
  // in other words: which of the format fields are set
  present.reset();
  formatChecked = 0;
  // first inser empty object 0, since object numbers start with 1
  vector<CatObject>::clear();
  CatObject s0 = {0,0,0,0,0,0,0,0,0,0};
  vector<CatObject>::push_back(s0);
  // open cat file
  ifstream catalog (catfile.c_str());
  if (catalog.fail()) {
    std::cerr << "Catalog: catalog file does not exist!" << endl;
    terminate();
  }
  catalog.clear();
  // read in cat file
  // 1) parse the format definition lines
  // 2) fill object information into SExCatObjects -> vector<CatObject>;
  string line;
  std::vector<std::string> column;
  while(getline(catalog, line)) {
    typedef boost::tokenizer<boost::char_separator<char> > Tok;
    // split entries at empty chars
    boost::char_separator<char> sep(" ");
    Tok tok(line, sep);
    // first of all we copy the token into string vector
    // though this is not too sophisticated it does not hurt and is more 
    // convenient
    column.clear();
    for(Tok::iterator tok_iter = tok.begin(); tok_iter != tok.end(); ++tok_iter)
      column.push_back(*tok_iter);
    // comment line: contains format definition at columns 2,3
    if (column[0].compare("#") == 0)
      setFormatField(column[2],atoi(column[1].c_str()));
    // at this point we should have a complete format definition: check it
    // from here on we expect the list of object to come along.
    else {
      if (!formatChecked) checkFormat();
      // then set up a true SExCatObject
      CatObject so;
      so.NUMBER = (unsigned int) atoi(column[format.NUMBER].c_str());
      // the sextractor corrdinates start with (1,1), ours with (0,0)
      so.XMIN = atoi(column[format.XMIN].c_str())-1;
      so.XMAX = atoi(column[format.XMAX].c_str())-1;
      so.YMIN = atoi(column[format.YMIN].c_str())-1;
      so.YMAX = atoi(column[format.YMAX].c_str())-1;
      so.XCENTROID = atof(column[format.XCENTROID].c_str())-1;
      so.YCENTROID = atof(column[format.YCENTROID].c_str())-1;
      so.FLUX = atof(column[format.FLUX].c_str());
      so.FLAGS = (unsigned char) atoi(column[format.FLAGS].c_str());
      if (present.test(9))
	so.CLASSIFIER = (data_t) atof(column[format.CLASSIFIER].c_str());
      else
	so.CLASSIFIER = (data_t) 0;
     // then push it on vector<CatObject>
      vector<CatObject>::push_back(so);
    }
  }
  catalog.close();
}

void Catalog::save(string catfile) const {
  // write header and the data section in SExtractor format
  ofstream catalog (catfile.c_str());
  catalog << "#  1 NUMBER" << endl;
  catalog << "#  2 XMIN_IMAGE" << endl;
  catalog << "#  3 XMAX_IMAGE" << endl;
  catalog << "#  4 YMIN_IMAGE" << endl;
  catalog << "#  5 YMAX_IMAGE" << endl;
  catalog << "#  6 XCENTROID_IMAGE" << endl;
  catalog << "#  7 YCENTROID_IMAGE" << endl;
  catalog << "#  8 FLUX" << endl;
  catalog << "#  9 FLAGS" << endl;
  if (present.test(9))
    catalog << "# 10 CLASSIFIER" << endl;
  // now all objects in Catalog
  for (unsigned int i=1; i<vector<CatObject>::size(); i++) {
    catalog << vector<CatObject>::operator[](i).NUMBER << " ";
    catalog << vector<CatObject>::operator[](i).XMIN + 1 << " ";
    catalog << vector<CatObject>::operator[](i).XMAX + 1 << " ";
    catalog << vector<CatObject>::operator[](i).YMIN + 1 << " ";
    catalog << vector<CatObject>::operator[](i).YMAX + 1 << " ";
    catalog << vector<CatObject>::operator[](i).XCENTROID + 1 << " ";
    catalog << vector<CatObject>::operator[](i).YCENTROID + 1 << " ";
    catalog << vector<CatObject>::operator[](i).FLUX << " ";
    catalog << (unsigned int) vector<CatObject>::operator[](i).FLAGS << " ";
    if (present.test(9))
      catalog << vector<CatObject>::operator[](i).CLASSIFIER << " ";
    catalog << endl;
  }
}

void Catalog::setFormatField(std::string type, unsigned short colnr) {
  if (type == "NUMBER" || type == "NR") {
    format.NUMBER = colnr - 1;
    present[0] = 1;
  }
  else if (type.find("XMIN") != string::npos) {
    format.XMIN = colnr - 1;
    present[1] = 1;
  }
  else if (type.find("XMAX") != string::npos) {
    format.XMAX = colnr - 1;
    present[2] = 1;
  }
  else if (type.find("YMIN") != string::npos) {
    format.YMIN = colnr - 1;
    present[3] = 1;
  }
  else if (type.find("YMAX") != string::npos) {
    format.YMAX = colnr - 1;
    present[4] = 1;
  }
  else if (type == "X" || type.find("X_") != string::npos || type.find("XWIN") != string::npos || type.find("XPEAK") != string::npos) {
    format.XCENTROID = colnr - 1;
    present[5] = 1;
  }
  else if (type == "Y" || type.find("Y_") != string::npos || type.find("YWIN") != string::npos || type.find("YPEAK") != string::npos) {
    format.YCENTROID = colnr - 1;
    present[6] = 1;
  } 
  else if (type.find("FLUX") != string::npos) {
    format.FLUX = colnr - 1;
    present[7] = 1;
  }
  else if (type == "FLAGS") {
    format.FLAGS = colnr - 1;
    present[8] = 1;
  }
  else if (type == "CLASS_STAR" || type == "CLASSIFIER") {
    format.CLASSIFIER = colnr - 1;
    present[9] = 1;
  }
}

bool Catalog::checkFormat() {
  // test if all but 9th bit are set
  // this means that CLASSIFIER is optional
  bitset<10> mask(0);
  mask[9] = 1;
  if ((present | mask).count() < 10) {
    cerr << "Catalog: mandatory catalog parameters are missing!" << endl;
    terminate();
  }
}
