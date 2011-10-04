#include "../../include/frame/Catalog.h"
#include "../../include/frame/Shapes.h"
#include <fstream>
#include <iostream>
#include <set>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <stdexcept>

using namespace shapelens;
using namespace std;

Catalog::Catalog() : map<unsigned long, CatObject>() {
  present.set();
}

Catalog::Catalog(string catfile) : map<unsigned long, CatObject>() {
  read(catfile);
}

void Catalog::read(string catfile) {
  // reset bitset that indicates which data are given in catalog
  // in other words: which of the format fields are set
  present.reset();
  formatChecked = true;
  
  // check if it's a fits table
  int status = 0;
  fitsfile* fptr = NULL;
  // yes, it is a FITS table
  try {
    fptr = IO::openFITSTable(catfile);
    // get column numbers for all required columns 
    setFormatFromFITSTable(fptr);
    checkFormat();
    // go thru all rows and grab the required columns for each object
    CatObject so;
    so.CLASSIFIER = so.PARENT = 0;
    long nrows = IO::getFITSTableRows(fptr);
    for (long i = 0; i < nrows; i++) {
      unsigned long id;
      IO::readFITSTableValue(fptr, i, format.ID, id);
      IO::readFITSTableValue(fptr, i, format.XMIN, so.XMIN);
      so.XMIN -= 1; // the sextractor coords start with (1,1), ours with (0,0)
      IO::readFITSTableValue(fptr, i, format.XMAX, so.XMAX);
      so.XMAX -= 1;
      IO::readFITSTableValue(fptr, i, format.YMIN, so.YMIN);
      so.YMIN -= 1;
      IO::readFITSTableValue(fptr, i, format.YMAX, so.YMAX);
      so.YMAX -= 1; 
      IO::readFITSTableValue(fptr, i, format.XCENTROID, so.XCENTROID);
      so.XCENTROID -= 1;
      IO::readFITSTableValue(fptr, i, format.YCENTROID, so.YCENTROID);
      so.YCENTROID -= 1;
      IO::readFITSTableValue(fptr, i, format.FLAGS, so.FLAGS);
      if (present.test(8))
	IO::readFITSTableValue(fptr, i, format.CLASSIFIER, so.CLASSIFIER);
      if (present.test(9))
	IO::readFITSTableValue(fptr, i, format.PARENT, so.PARENT);
      // then store it in map
      Catalog::insert(make_pair(id,so));
    }
    IO::closeFITSFile(fptr);
  }

  // no, it's an ASCII file
  catch (std::runtime_error) {
    // open cat file
    ifstream catalog (catfile.c_str());
    if (catalog.fail())
      throw std::invalid_argument("Catalog: catalog file does not exist!");
    catalog.clear();

    // read in cat file
    // 1) parse the format definition lines
    // 2) fill object information into SExCatObjects -> map<unsigned long, CatObject>;
    string line;
    std::vector<std::string> column;
    while(getline(catalog, line)) {
      typedef boost::tokenizer<boost::char_separator<char> > Tok;
      // split entries at empty chars
      boost::char_separator<char> sep(" \t");
      Tok tok(line, sep);
      // first of all we copy the token into string vector
      // though this is not too sophisticated it does not hurt and is more 
      // convenient
      column.clear();
      for(Tok::iterator tok_iter = tok.begin(); tok_iter != tok.end(); ++tok_iter)
	column.push_back(*tok_iter);
      
      // exclude empty lines
      if (column.size() > 2) { // at least # COLUMN NAME needs to be present 
	// comment line: contains format definition at columns 2,3
	if (column[0].compare("#") == 0)
	  setFormatField(column[2],boost::lexical_cast<unsigned short>(column[1].c_str()));
	// at this point we should have a complete format definition: check it
	// from here on we expect the list of object to come along.
	else {
	  if (!formatChecked) checkFormat();
	  if (column.size() < present.count())
	    throw std::runtime_error("Catalog: not all necessary column provided in line:\n" + line);
	  // then set up a true SExCatObject
	  CatObject so;
	  so.CLASSIFIER = so.PARENT = 0;
	  unsigned long id = boost::lexical_cast<unsigned long>(column[format.ID].c_str());
	  // the sextractor corrdinates start with (1,1), ours with (0,0)
	  so.XMIN = boost::lexical_cast<int>(column[format.XMIN].c_str())-1;
	  so.XMAX = boost::lexical_cast<int>(column[format.XMAX].c_str())-1;
	  so.YMIN = boost::lexical_cast<int>(column[format.YMIN].c_str())-1;
	  so.YMAX = boost::lexical_cast<int>(column[format.YMAX].c_str())-1;
	  // as Catalog pixels start with (1,1) and ours with (0,0), 
	  // we need to subtract 1
	  so.XCENTROID = boost::lexical_cast<data_t>(column[format.XCENTROID].c_str())-1;
	  so.YCENTROID = boost::lexical_cast<data_t>(column[format.YCENTROID].c_str())-1;
	  so.FLAGS = (unsigned char) boost::lexical_cast<unsigned int>(column[format.FLAGS].c_str());
	  if (present.test(8))
	    so.CLASSIFIER = boost::lexical_cast<data_t>(column[format.CLASSIFIER].c_str());
	  if (present.test(9))
	    so.PARENT = boost::lexical_cast<unsigned long>(column[format.PARENT].c_str());
	  // then store it in map
	  Catalog::insert(make_pair(id,so));
	}
      }
    }
    catalog.close();
  }
}

void Catalog::save(string catfile) const {
  // write header and the data section in SExtractor format
  ofstream catalog (catfile.c_str());
  catalog << "#  1 ID" << endl;
  catalog << "#  2 XMIN_IMAGE" << endl;
  catalog << "#  3 XMAX_IMAGE" << endl;
  catalog << "#  4 YMIN_IMAGE" << endl;
  catalog << "#  5 YMAX_IMAGE" << endl;
  catalog << "#  6 X_IMAGE" << endl;
  catalog << "#  7 Y_IMAGE" << endl;
  catalog << "#  8 FLAGS" << endl;
  if (present.test(8))
    catalog << "# 9 CLASSIFIER" << endl;
  if (present.test(9)) {
    if (present.test(8))
      catalog << "# 10 PARENT" << endl;
    else
      catalog << "# 9 PARENT" << endl;
  }
  // now all objects in Catalog
  Catalog::const_iterator iter;
  for (iter = Catalog::begin(); iter != Catalog::end(); iter++) {
    catalog << iter->first << " ";
    catalog << iter->second.XMIN + 1 << " ";
    catalog << iter->second.XMAX + 1 << " ";
    catalog << iter->second.YMIN + 1 << " ";
    catalog << iter->second.YMAX + 1 << " ";
    catalog << iter->second.XCENTROID + 1 << " ";
    catalog << iter->second.YCENTROID + 1 << " ";
    catalog << (unsigned int) iter->second.FLAGS << " ";
    if (present.test(8))
      catalog << iter->second.CLASSIFIER << " ";
    if (present.test(9))
      catalog << iter->second.PARENT << " ";
    catalog << endl;
  }
}

void Catalog::setFormatField(std::string type, unsigned short colnr) {
  if (present[0] == 0 && (type == "ID" || type == "NUMBER" || type == "NR")) {
    format.ID = colnr - 1;
    present[0] = 1;
  }
  else if (present[1] == 0 && type.find("XMIN") == 0) {
    format.XMIN = colnr - 1;
    present[1] = 1;
  }
  else if (present[2] == 0 && type.find("XMAX") == 0) {
    format.XMAX = colnr - 1;
    present[2] = 1;
  }
  else if (present[3] == 0 && type.find("YMIN") == 0) {
    format.YMIN = colnr - 1;
    present[3] = 1;
  }
  else if (present[4] == 0 && type.find("YMAX") == 0) {
    format.YMAX = colnr - 1;
    present[4] = 1;
  }
  else if (present[5] == 0 && (type == "X" || type.find("X_") == 0 || type.find("XWIN") == 0 || type.find("XPEAK") == 0 || type.find("XCENTROID") == 0)) {
    format.XCENTROID = colnr - 1;
    present[5] = 1;
  }
  else if (present[6] == 0 && (type == "Y" || type.find("Y_") == 0 || type.find("YWIN") == 0 || type.find("YPEAK") == 0 || type.find("YCENTROID") == 0)) {
    format.YCENTROID = colnr - 1;
    present[6] = 1;
  } 
  else if (present[7] == 0 && type == "FLAGS") {
    format.FLAGS = colnr - 1;
    present[7] = 1;
  }
  else if (present[8] == 0 && (type == "CLASS_STAR" || type == "CLASS" || type == "CLASSIFIER")) {
    format.CLASSIFIER = colnr - 1;
    present[8] = 1;
  }
  else if (present[9] == 0 && (type == "PARENT" || type == "VECTOR_ASSOC")) {
    format.PARENT = colnr - 1;
    present[9] = 1;
  }
}

void Catalog::setFormatFromFITSTable(fitsfile* fptr) {
  try {
    format.ID = IO::getFITSTableColumnNumber(fptr, "ID");
    present[0] = 1;
  } catch (std::invalid_argument) {
    try {
      format.ID = IO::getFITSTableColumnNumber(fptr, "NUMBER");
      present[0] = 1;
    } catch (std::invalid_argument) {
      try {
	format.ID = IO::getFITSTableColumnNumber(fptr, "NR");
	present[0] = 1;
      } catch (std::invalid_argument) {}
    }
  }
  try {
    format.XMIN = IO::getFITSTableColumnNumber(fptr, "XMIN");
    present[1] = 1;
  } catch (std::invalid_argument) {
    try {
      format.XMIN = IO::getFITSTableColumnNumber(fptr, "XMIN_IMAGE");
      present[1] = 1;
    } catch  (std::invalid_argument) {}
  }
  try {
    format.XMAX = IO::getFITSTableColumnNumber(fptr, "XMAX");
    present[2] = 1;
  } catch (std::invalid_argument) {
    try {
      format.XMAX = IO::getFITSTableColumnNumber(fptr, "XMAX_IMAGE");
      present[2] = 1;
    } catch (std::invalid_argument) {}
  }
  try {
    format.YMIN = IO::getFITSTableColumnNumber(fptr, "YMIN");
    present[3] = 1;
  } catch (std::invalid_argument) {
    try {
      format.YMIN = IO::getFITSTableColumnNumber(fptr, "YMIN_IMAGE");
      present[3] = 1;
    } catch (std::invalid_argument) {}
  }
  try {
    format.YMAX = IO::getFITSTableColumnNumber(fptr, "YMAX");
    present[4] = 1;
  } catch (std::invalid_argument) {
    try {
      format.YMAX = IO::getFITSTableColumnNumber(fptr, "YMAX_IMAGE");
      present[4] = 1;
    } catch (std::invalid_argument) {}
  }
  try {
    format.XCENTROID = IO::getFITSTableColumnNumber(fptr, "X");
    present[5] = 1;
  } catch (std::invalid_argument) {
    try {
      format.XCENTROID = IO::getFITSTableColumnNumber(fptr, "X_IMAGE");
      present[5] = 1;
    } catch (std::invalid_argument) {
      try {
	format.XCENTROID = IO::getFITSTableColumnNumber(fptr, "XWIN_IMAGE");
	present[5] = 1;
      } catch (std::invalid_argument) {
	try {
	  format.XCENTROID = IO::getFITSTableColumnNumber(fptr, "XCENTROID");
	  present[5] = 1;
	} catch (std::invalid_argument) {
	  try {
	    format.XCENTROID = IO::getFITSTableColumnNumber(fptr, "XPEAK");
	    present[5] = 1;
	  } catch (std::invalid_argument) {}}
      }
    }
  }
  try {
    format.YCENTROID = IO::getFITSTableColumnNumber(fptr, "Y");
    present[6] = 1;
  } catch (std::invalid_argument) {
    try {
      format.YCENTROID = IO::getFITSTableColumnNumber(fptr, "Y_IMAGE");
      present[6] = 1;
    } catch (std::invalid_argument) {
      try {
	format.YCENTROID = IO::getFITSTableColumnNumber(fptr, "YWIN_IMAGE");
	present[6] = 1;
      } catch (std::invalid_argument) {
	try {
	  format.YCENTROID = IO::getFITSTableColumnNumber(fptr, "YCENTROID");
	  present[6] = 1;
	} catch (std::invalid_argument) {
	  try {
	    format.YCENTROID = IO::getFITSTableColumnNumber(fptr, "YPEAK");
	    present[6] = 1;
	  } catch (std::invalid_argument) {}
	}
      }
    }
  }
  try {
    format.FLAGS = IO::getFITSTableColumnNumber(fptr, "FLAGS");
    present[7] = 1;
  } catch (std::invalid_argument) {}
  try {
    format.CLASSIFIER = IO::getFITSTableColumnNumber(fptr, "CLASSIFIER");
    present[8] = 1;
  } catch (std::invalid_argument) {
    try {
    format.CLASSIFIER = IO::getFITSTableColumnNumber(fptr, "CLASS_STAR");
    present[8] = 1;
    } catch (std::invalid_argument) {}
  }
  try {
    format.PARENT = IO::getFITSTableColumnNumber(fptr, "PARENT");
    present[9] = 1;
  } catch (std::invalid_argument) {
    try {
      format.PARENT = IO::getFITSTableColumnNumber(fptr, "VECTOR_ASSOC");
      present[9] = 1;
    } catch (std::invalid_argument) {}
  }
}

bool Catalog::checkFormat() {
  // test if all but 8th and 9th bit are set
  // this means that CLASSIFIER and PARENT are optional
  bitset<10> mask(0);
  mask[8] = mask[9] = 1;
  if ((present | mask).count() < 10) {
    cerr << "Catalog: mandatory catalog parameters are missing!" << endl;
    terminate();
  }
  else
    formatChecked = 1;
}
  

int Catalog::round(data_t x) {
  if (x-floor(x)>=0.5)
    return (int) ceil(x);
  else
    return (int) floor(x);
}

void Catalog::apply(const CoordinateTransformation& C) {
  Point<data_t> P;
  Rectangle<data_t> bb;
  for (Catalog::iterator iter = Catalog::begin(); iter != Catalog::end(); iter++) {
    CatObject& ci = iter->second;
    // transform centroid
    P(0) = ci.XCENTROID;
    P(1) = ci.YCENTROID;
    C.transform(P);
    ci.XCENTROID = P(0);
    ci.YCENTROID = P(1);
    // transform bounding box of (XMIN/YMIN) .. (XMAX/YMAX)
    bb.ll(0) = ci.XMIN;
    bb.ll(1) = ci.YMIN;
    bb.tr(0) = ci.XMAX;
    bb.tr(1) = ci.YMAX;
    bb.apply(C);
    ci.XMIN = (int) round(bb.ll(0));
    ci.YMIN = (int) round(bb.ll(1));
    ci.XMAX = (int) round(bb.tr(0));
    ci.YMAX = (int) round(bb.tr(1));
  }
}
  
Catalog Catalog::operator+ (const Catalog& c) {
  Catalog newCat = *this;
  for (Catalog::const_iterator iter = c.begin(); iter != c.end(); iter++)
    newCat[iter->first] = iter->second;
}

void Catalog::operator+= (const Catalog& c) {
  for (Catalog::const_iterator iter = c.begin(); iter != c.end(); iter++)
    Catalog::operator[](iter->first) = iter->second;
}

Catalog Catalog::operator- (const Catalog& c) {
  Catalog newCat = *this;
  for (Catalog::const_iterator iter = c.begin(); iter != c.end(); iter++) {
    Catalog::iterator key = newCat.find(iter->first);
    if (key != newCat.end())
      newCat.erase(key);
  }
}

void Catalog::operator-= (const Catalog& c) {
  for (Catalog::const_iterator iter = c.begin(); iter != c.end(); iter++) {
    Catalog::iterator key = Catalog::find(iter->first);
    if (key != Catalog::end())
      Catalog::erase(key);
  }
}

