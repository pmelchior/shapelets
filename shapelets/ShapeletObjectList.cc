#include <ShapeletObjectList.h>
#include <fstream>

using namespace std;

bool alwaysTrue(ShapeletObject& so) {
  return 1;
}

ShapeletObjectList::ShapeletObjectList(string filename) : list<ShapeletObject*>() {
  readListFile(filename,&alwaysTrue);
}

ShapeletObjectList::ShapeletObjectList(string filename, bool (* selectionFunction) (ShapeletObject&)) : list<ShapeletObject*>() {
  readListFile(filename,selectionFunction);
}

ShapeletObjectList::~ShapeletObjectList() {
  // we have to destroy the ShapeletObject via its pointer
  list<ShapeletObject*>::iterator iter;
  for (iter = list<ShapeletObject*>::begin(); iter != list<ShapeletObject*>::end(); iter++) 
    delete *iter;
}
void ShapeletObjectList::readListFile(string filename, bool (* selectionFunction) (ShapeletObject&)) {
  // open file with list of ShapeletObjects
  ifstream listfile (filename.c_str());
  if (listfile.fail()) {
    cout << "ShapeletObjectList: list file does not exists!" << endl;
    terminate();
  }
  // read in list file
  string sifname;
  while(getline(listfile, sifname)) {
    ShapeletObject* sobj = new ShapeletObject(sifname);
    // if selectionFunction return 1 for this sobj,
    // append it to list
    if ((*selectionFunction)(*sobj))
      list<ShapeletObject*>::push_back(sobj);
    else
      delete sobj;
  }
}
