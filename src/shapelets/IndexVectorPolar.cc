#include "../../include/shapelets/IndexVectorPolar.h"

using namespace shapelens;
using namespace std;

IndexVectorPolar::IndexVectorPolar () :
  nmax(0), pMatrix(1,1) {
  computeIndexMaps();
}

IndexVectorPolar::IndexVectorPolar (unsigned int innmax) :
  nmax(innmax), pMatrix(innmax+1, innmax+1) {
  computeIndexMaps();
}

void IndexVectorPolar::setNMax(unsigned int innmax) {
  if (innmax != nmax) {
    nmax = innmax;
    polar.clear();
    pMatrix.resize(nmax+1,nmax+1);
    pMatrix.clear();
    computeIndexMaps();
  }
}  

int IndexVectorPolar::getNMax() const {
  return nmax;
}

int IndexVectorPolar::getNCoeffs() const {
  return (nmax+1)*(nmax+2)/2;
}

int IndexVectorPolar::mIndex(unsigned int n, int m) const {
  return (m+n)/2;
}

void IndexVectorPolar::computeIndexMaps() {
  unsigned int i = 0;
  for (int n=0; n<= nmax; n++) {
    for (int n1i=0; n1i <=n; n1i++) {
      polar[i] = pair<unsigned int, int>(n,2*n1i - n);
      pMatrix(n,mIndex(n,2*n1i - n)) = i;
      i++;
    }
  }
}

int IndexVectorPolar::getState1(unsigned int i) const {
  std::map<unsigned int, std::pair<unsigned int, int> >::const_iterator pIter = polar.find(i);
  if (pIter != polar.end())
    return pIter->second.first;
  else {
    std::cerr << "IndexVector: index " << i << " does not exist" << std::endl;
    std::terminate();
  }
}

int IndexVectorPolar::getState2(unsigned int i) const {
  std::map<unsigned int, std::pair<unsigned int, int> >::const_iterator pIter = polar.find(i);
  if (pIter != polar.end())
    return pIter->second.second;
  else {
    std::cerr << "IndexVector: index " << i << " does not exist" << std::endl;
    std::terminate();
  }
}

int IndexVectorPolar::getIndex1(unsigned int i) const {
  return IndexVectorPolar::getState1(i);
}

int IndexVectorPolar::getIndex2(unsigned int i) const {
  std::map<unsigned int, std::pair<unsigned int, int> >::const_iterator pIter = polar.find(i);
  if (pIter != polar.end())
    return mIndex(pIter->second.first,pIter->second.second);
  else {
    std::cerr << "IndexVector: index " << i << " does not exist" << std::endl;
    std::terminate();
  }
}

int IndexVectorPolar::getIndex(int n, int m) const {
  if (n <= nmax && abs(m) <= n)
    return pMatrix(n,mIndex(n,m));
  else {
    std::cerr << "IndexVector: eigenstate combination " << n << "/" << m << " illegal!" << std::endl;
    std::terminate();
  }
}

