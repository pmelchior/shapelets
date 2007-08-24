#include <ShapeletObjectList.h>
#include <fstream>

using namespace std;

bool alwaysTrue(ShapeletObject& so) {
  return 1;
}

ShapeletObjectList::ShapeletObjectList(string filename) : vector<boost::shared_ptr<ShapeletObject> >() {
  readListFile(filename,&alwaysTrue);
}

ShapeletObjectList::ShapeletObjectList(string filename, bool (* selectionFunction) (ShapeletObject&)) : vector<boost::shared_ptr<ShapeletObject> >() {
  readListFile(filename,selectionFunction);
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
    boost::shared_ptr<ShapeletObject> safePointer (new ShapeletObject(sifname));
    // if selectionFunction return 1 for this sobj,
    // append it to list
    if ((*selectionFunction)(*safePointer))
      vector<boost::shared_ptr<ShapeletObject> >::push_back(safePointer);
  }
}
