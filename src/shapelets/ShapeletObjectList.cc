#include <shapelets/ShapeletObjectList.h>
#include <fstream>
#include <list>
#include <cmath>

using namespace shapelens;
using namespace std;

ShapeletObjectList::ShapeletObjectList() : vector<boost::shared_ptr<ShapeletObject> >() {
}

ShapeletObjectList::ShapeletObjectList(string filename) : vector<boost::shared_ptr<ShapeletObject> >() {
  readListFile(filename,NULL,NULL);
}

ShapeletObjectList::ShapeletObjectList(string filename, bool (* selectionFunction) (ShapeletObject&, void*), void* p) : vector<boost::shared_ptr<ShapeletObject> >() {
  readListFile(filename,selectionFunction,p);
}

void ShapeletObjectList::readListFile(string filename, bool (* selectionFunction) (ShapeletObject&, void*), void* p) {
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
    // if there is a selectionFunction and it returns 1 for this sobj,
    // append it to list
    if (selectionFunction != NULL) {
      if ((*selectionFunction)(*safePointer,p))
	vector<boost::shared_ptr<ShapeletObject> >::push_back(safePointer);
    }else // otherwise just append it
      vector<boost::shared_ptr<ShapeletObject> >::push_back(safePointer);
  }
}

void ShapeletObjectList::average(CoefficientVector<data_t>& mean, CoefficientVector<data_t>& std_mean, data_t& beta) {
  average(mean,std_mean,beta,NULL);
}

void ShapeletObjectList::average(CoefficientVector<data_t>& mean, CoefficientVector<data_t>& std_mean, data_t& beta, data_t (* weightFunction) (ShapeletObject&, void*), void* p) {
  // set up two empty matrices for mean and std_mean
  // as they are easier to resize
  NumMatrix<data_t> meanMatrix(1,1), stdMatrix(1,1);
  beta = 0;
  data_t sum_weights = 0, sum_weights2 = 0;
  int nmax=0, n1, n2;

  // go through all ShapeletObjects
  for (ShapeletObjectList::iterator iter = ShapeletObjectList::begin(); iter != ShapeletObjectList::end(); iter++) {
    const CoefficientVector<data_t>& coeffs = (*iter)->getCoeffs();
    const IndexVector& nVector = coeffs.getIndexVector();
    data_t weight = 1;
    if (weightFunction != NULL)
      weight = (*weightFunction)(*(*iter), p);
    // if new coeff matrix is bigger than current matrices:
    // expand them
    if (coeffs.getNMax() > nmax) {
      nmax = coeffs.getNMax();
      meanMatrix.resize_clear(nmax+1, nmax+1);
      stdMatrix.resize_clear(nmax+1,nmax+1);
    }
    // go thru all coeffs
    for (unsigned int i=0; i < nVector.getNCoeffs(); i++) {
      n1 = nVector.getState1(i);
      n2 = nVector.getState2(i);
      meanMatrix(n1,n2) += weight*coeffs(i);
      stdMatrix(n1,n2) += weight*coeffs(i)*coeffs(i);
    }
    beta += weight*((*iter)->getBeta());
    sum_weights += weight;
    sum_weights2 += weight*weight;
  }
  // compute average of beta and all coeffs
  beta /= sum_weights;
  
  // set up mean and std_mean with appropriate nmax
  mean.setNMax(nmax);
  std_mean.setNMax(nmax);
  const IndexVector& nVector = mean.getIndexVector();
  for (unsigned int i=0; i < nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    std_mean(i) = sqrt((stdMatrix(n1,n2)*sum_weights - meanMatrix(n1,n2)*meanMatrix(n1,n2)) / (sum_weights*sum_weights - sum_weights2));
    mean(i) = meanMatrix(n1,n2) / sum_weights;
  }
}

ShapeletObjectList ShapeletObjectList::select(bool (* selectionFunction) (ShapeletObject&, void*), void* p) {
  ShapeletObjectList selection;
  for(vector<boost::shared_ptr<ShapeletObject> >::iterator iter = vector<boost::shared_ptr<ShapeletObject> >::begin(); iter != vector<boost::shared_ptr<ShapeletObject> >::end(); iter++)
    if ((*selectionFunction)(*(*iter),p))
      selection.push_back(*iter);				
  return selection;
}		
