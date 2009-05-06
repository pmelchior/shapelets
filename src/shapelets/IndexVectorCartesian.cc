#include "../../include/shapelets/IndexVectorCartesian.h"

using namespace shapelens;
using namespace std;

IndexVectorCartesian::IndexVectorCartesian () :
  nmax(0), cMatrix(1,1) {
  computeIndexMaps();
}

IndexVectorCartesian::IndexVectorCartesian (unsigned int innmax) :
  nmax(innmax), cMatrix(innmax+1, innmax+1) {
  computeIndexMaps();
}

void IndexVectorCartesian::setNMax(unsigned int innmax) {
  if (innmax != nmax) {
    nmax = innmax;
    cartesian.clear();
    cMatrix.resize(nmax+1,nmax+1);
    cMatrix.clear();
    computeIndexMaps();
  }
}  

int IndexVectorCartesian::getNMax() const {
  return nmax;
}

int IndexVectorCartesian::getNCoeffs() const {
  return (nmax+1)*(nmax+2)/2;
}

void IndexVectorCartesian::computeIndexMaps() {
  unsigned int i = 0;
  for (int n=0; n<= nmax; n++) {
    for (int n1i=0; n1i <=n; n1i++) {
      cartesian[i] = pair<unsigned int, unsigned int>(n1i,n-n1i);
      cMatrix(n1i, n-n1i) = i;
      i++;
    }
  }
}

int IndexVectorCartesian::getState1(unsigned int i) const {
  std::map<unsigned int, std::pair<unsigned int, unsigned int> >::const_iterator cIter = cartesian.find(i);
  if (cIter != cartesian.end())
    return cIter->second.first;
  else {
    std::cerr << "IndexVector: index " << i << " does not exist" << std::endl;
    std::terminate();
  }
}

int IndexVectorCartesian::getState2(unsigned int i) const {
  std::map<unsigned int, std::pair<unsigned int, unsigned int> >::const_iterator cIter = cartesian.find(i);
  if (cIter != cartesian.end())
    return cIter->second.second;
  else {
    std::cerr << "IndexVector: index " << i << " does not exist" << std::endl;
    std::terminate();
  }
}

int IndexVectorCartesian::getIndex1(unsigned int i) const {
  return IndexVectorCartesian::getState1(i);
}

int IndexVectorCartesian::getIndex2(unsigned int i) const {
  return IndexVectorCartesian::getState2(i);
}

int IndexVectorCartesian::getIndex(int n1, int n2) const {
  if (n1 + n2 <= nmax)
    return cMatrix(n1,n2);
  else {
    std::cerr << "IndexVector: eigenstate combination " << n1 << "/" << n2 << " illegal!" << std::endl;
    std::terminate();
  }
}
