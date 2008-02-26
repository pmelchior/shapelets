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

void ShapeletObjectList::average(CoefficientVector<data_t>& mean, CoefficientVector<data_t>& std_mean, data_t& beta) {
  average(mean,std_mean,beta,&alwaysOne);
}

void ShapeletObjectList::average(CoefficientVector<data_t>& mean, CoefficientVector<data_t>& std_mean, data_t& beta, data_t (* weightFunction) (ShapeletObject&)) {
  // set up two empty matrices for mean and std_mean
  // as they are easier to resize
  NumMatrix<data_t> meanMatrix, stdMatrix;
  beta = 0;
  data_t sum_weights = 0, sum_weights2 = 0;
  unsigned int nmax, n1, n2;

  // go through all ShapeletObjects
  for (ShapeletObjectList::iterator iter = ShapeletObjectList::begin(); iter != ShapeletObjectList::end(); iter++) {
    const CoefficientVector<data_t>& coeffs = (*iter)->getCoeffs();
    const IndexVector& nVector = coeffs.getIndexVector();
    data_t weight = (*weightFunction)(*(*iter));
    // if new coeff matrix is bigger than current matrices:
    // expand them
    if (coeffs.getNMax() > nmax) {
      nmax = coeffs.getNMax();
      meanMatrix.resize_clear(nmax+1, nmax+1);
      stdMatrix.resize_clear(nmax+1,nmax+1);
    }
    // go thru all coeffs
    for (unsigned int i=0; i < nVector.getNCoeffs(); i++) {
      n1 = nVector.getN1(i);
      n2 = nVector.getN2(i);
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
    n1 = nVector.getN1(i);
    n2 = nVector.getN2(i);
    std_mean(i) = sqrt((stdMatrix(n1,n2)*sum_weights - meanMatrix(n1,n2)*meanMatrix(n1,n2)) / (sum_weights*sum_weights - sum_weights2));
    mean(i) = meanMatrix(n1,n2) / sum_weights;
  }
}

bool ShapeletObjectList::checkSIFFile(ShapeletObject& so, std::string siffile) {
  // check coeffs and beta for correctness
  const CoefficientVector<data_t>& coeffs = so.getCoeffs();
  const IndexVector& nVector = coeffs.getIndexVector();
  bool ok = 1;
  for (unsigned int i=0; i < nVector.getNCoeffs(); i++) {
    if (isinf(coeffs(i)) || isnan(coeffs(i)) || fabs(coeffs(i))> 1e20) {
      ok = 0;
      break;
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
