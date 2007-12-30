#include <shapelets/ShapeletObjectList.h>
#include <fstream>
#include <list>
#include <cmath>

using namespace std;

bool alwaysTrue(ShapeletObject& so) {
  return 1;
}

data_t alwaysOne(ShapeletObject& so) {
  return 1;
}

ShapeletObjectList::ShapeletObjectList() : vector<boost::shared_ptr<ShapeletObject> >() {
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
    if ((*selectionFunction)(*safePointer) && checkSIFFile(*safePointer,sifname))
      vector<boost::shared_ptr<ShapeletObject> >::push_back(safePointer);
  }
}

void ShapeletObjectList::average(NumMatrix<data_t>& mean, NumMatrix<data_t>& std_mean, data_t& beta) {
  average(mean,std_mean,beta,&alwaysOne);
}

void ShapeletObjectList::average(NumMatrix<data_t>& mean, NumMatrix<data_t>& std_mean, data_t& beta, data_t (* weightFunction) (ShapeletObject&)) {
  //clear average and std_mean
  mean.clear();
  std_mean.clear();
  beta = 0;
  
  data_t sum_weights = 0, sum_weights2 = 0;

  // go through all ShapeletObjects
  for (ShapeletObjectList::iterator iter = ShapeletObjectList::begin(); iter != ShapeletObjectList::end(); iter++) {
    const NumMatrix<data_t>& coeffs = (*iter)->getCoeffs();
    data_t weight = (*weightFunction)(*(*iter));
    // if new coeff matrix is bigger than current average matrix
    // expand average
    if (mean.getRows() < coeffs.getRows() || mean.getColumns() < coeffs.getColumns()) {
      mean.resize_clear(coeffs.getRows(),coeffs.getColumns());
      std_mean.resize_clear(coeffs.getRows(),coeffs.getColumns());
    }
    for (int i=0; i<coeffs.getRows(); i++) {
      for (int j=0; j<coeffs.getColumns()-i; j++) {
	mean(i,j) += weight*coeffs(i,j);
	std_mean(i,j) += weight*coeffs(i,j)*coeffs(i,j);
      }
    }
    beta += weight*((*iter)->getBeta());
    sum_weights += weight;
    sum_weights2 += weight*weight;
  }
  // compute average of beta and all coeffs
  beta /= sum_weights;
  
  for (int i=0; i<mean.getRows(); i++) {
    for (int j=0; j<mean.getColumns()-i; j++) {
      std_mean(i,j) = sqrt((std_mean(i,j)*sum_weights - mean(i,j)*mean(i,j)) /
			   (sum_weights*sum_weights - sum_weights2));
      mean(i,j) /= sum_weights;
    }
  }
}

bool ShapeletObjectList::checkSIFFile(ShapeletObject& so, std::string siffile) {
  // check coeffs and beta for correctness
  const NumMatrix<data_t>& coeffs = so.getCoeffs();
  bool ok = 1;
  for (int i=0; i<coeffs.getRows(); i++) {
    for (int j=0; j<coeffs.getColumns(); j++) {
      if (isinf(coeffs(i,j)) || isnan(coeffs(i,j)) || fabs(coeffs(i,j))> 1e20) {
	ok = 0;
	break;
      }
    }
  }
  if (!ok)
    std::cerr << siffile << " contains errors, thus rejected from list." << std::cout;
  return ok;
}

ShapeletObjectList ShapeletObjectList::select(bool (* selectionFunction) (ShapeletObject&)) {
  ShapeletObjectList selection;
  for(vector<boost::shared_ptr<ShapeletObject> >::iterator iter = vector<boost::shared_ptr<ShapeletObject> >::begin(); iter != vector<boost::shared_ptr<ShapeletObject> >::end(); iter++)
    if ((*selectionFunction)(*(*iter)))
      selection.push_back(*iter);				
  return selection;
}		
