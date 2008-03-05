#include <frame/Catalog.h>
#include <frame/Point2D.h>
#include <fstream>
#include <iostream>
#include <set>
#include <boost/tokenizer.hpp>
#include <gsl/gsl_math.h>

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
  formatChecked = 0;
  // open cat file
  ifstream catalog (catfile.c_str());
  if (catalog.fail()) {
    std::cerr << "Catalog: catalog file does not exist!" << endl;
    terminate();
  }
  catalog.clear();
  // read in cat file
  // 1) parse the format definition lines
  // 2) fill object information into SExCatObjects -> map<unsigned long, CatObject>;
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
      unsigned long id = (unsigned long) atoi(column[format.ID].c_str());
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
      if (present.test(10))
	so.PARENT = (unsigned long) atoi(column[format.PARENT].c_str());
      else
	so.PARENT = (unsigned long)0;
      // then store it in map
      Catalog::operator[](id) = so;
    }
  }
  catalog.close();
}

void Catalog::save(string catfile) const {
  // write header and the data section in SExtractor format
  ofstream catalog (catfile.c_str());
  catalog << "#  1 ID" << endl;
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
  if (present.test(10)) {
    if (present.test(9))
      catalog << "# 11 PARENT" << endl;
    else
      catalog << "# 10 PARENT" << endl;
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
    catalog << iter->second.FLUX << " ";
    catalog << (unsigned int) iter->second.FLAGS << " ";
    if (present.test(9))
      catalog << iter->second.CLASSIFIER << " ";
    if (present.test(10))
      catalog << iter->second.PARENT << " ";
    catalog << endl;
  }
}

void Catalog::setFormatField(std::string type, unsigned short colnr) {
  if (type == "ID" || type == "NUMBER" || type == "NR") {
    format.ID = colnr - 1;
    present[0] = 1;
  }
  else if (type.find("XMIN") == 0) {
    format.XMIN = colnr - 1;
    present[1] = 1;
  }
  else if (type.find("XMAX") == 0) {
    format.XMAX = colnr - 1;
    present[2] = 1;
  }
  else if (type.find("YMIN") == 0) {
    format.YMIN = colnr - 1;
    present[3] = 1;
  }
  else if (type.find("YMAX") == 0) {
    format.YMAX = colnr - 1;
    present[4] = 1;
  }
  else if (type == "X" || type.find("X_") == 0 || type.find("XWIN") == 0 || type.find("XPEAK") == 0 || type.find("XCENTROID") == 0) {
    format.XCENTROID = colnr - 1;
    present[5] = 1;
  }
  else if (type == "Y" || type.find("Y_") == 0 || type.find("YWIN") == 0 || type.find("YPEAK") == 0 || type.find("YCENTROID") == 0) {
    format.YCENTROID = colnr - 1;
    present[6] = 1;
  } 
  else if (type.find("FLUX") == 0) {
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
  else if (type == "PARENT" || type == "VECTOR_ASSOC") {
    format.PARENT = colnr - 1;
    present[10] = 1;
  }
}

bool Catalog::checkFormat() {
  // test if all but 9th bit are set
  // this means that CLASSIFIER is optional
  bitset<11> mask(0);
  mask[9] = mask[10] = 1;
  if ((present | mask).count() < 11) {
    cerr << "Catalog: mandatory catalog parameters are missing!" << endl;
    terminate();
  }
  else
    formatChecked = 1;
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

// struct matchProps {
//   Point2D centroid;
//   data_t flux;
//   data_t offset;
// };

// struct matchBundle {
//   const Point2D& refCentroid;
//   const Catalog& cat;
//   map<ulong, matchProps>& matches;
// };

// bool insertCatMatch(ulong id, void* params) {
//   matchBundle* bundle = (matchBundle*) params;
//   Catalog::const_iterator iter = bundle->cat.find(id);
//   const Point2D& ref = bundle->refCentroid;
//   matchProps props = { Point2D(iter->second.XCENTROID, iter->second.YCENTROID), iter->second.FLUX, 0};
//   props.offset = sqrt(gsl_pow_2(props.centroid(0) - ref(0)) + gsl_pow_2(props.centroid(1) - ref(1)));
//   bundle->matches[id] = props;
//   return true;
// }

// Catalog Catalog::operator* (const Catalog& c) {
//   RTree<unsigned long,unsigned long,2,data_t> thisTree, cTree;
//   buildRTree(thisTree,*this);
//   buildRTree(cTree,c);
  
//   // set up structures
//   multimap<ulong, ulong> mapping;
//   map<ulong, matchProps> matchesC, matchesThis;
//   int foundC, foundThis;
//   Point2D refCentroidC, refCentroidThis;
//   matchBundle mbC { refCentroidThis, c, matchesC};
//   matchBundle mbThis { refCentroidThis, *this, matchesThis};

//   // iterate through objects of the smaller catalog
//   Catalog::Rectangle<ulong> searchRectC, searchRectThis;
//   if (Catalog::size() < c.size()) {
//     for (Catalog::const_iterator iter = Catalog::begin(); iter != Catalog::end(); iter++) {
//       refCentroidThis(0) = iter->second.XCENTROID;
//       refCentroidThis(1) = iter->second.YCENTROID;
//       matchesC.clear();
//       searchRectThis.setCoords(iter);
//       foundC = cTree.Search(searchRectThis.getMin(), searchRectThis.getMax(), &insertCatMatch, &matchesC);
//     }
//   } else {
//     for (Catalog::const_iterator iter = c.begin(); iter != c.end(); iter++) {
//     }
//   }
// }

