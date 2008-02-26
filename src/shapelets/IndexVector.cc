#include <shapelets/IndexVector.h>
#include <algorithm>

using namespace std;

IndexVector::IndexVector () {
  nmax = 0;
}

IndexVector::IndexVector (int nmax) {
  setNMax(nmax);
}

void IndexVector::setNMax(int innmax) {
  if (innmax != nmax) {
    nmax = innmax;
    cartesian.clear();
    polar.clear();
    cMatrix.resize(nmax+1,nmax+1);
    cMatrix.clear();
    pMatrix.resize(nmax+1,nmax+1);
    pMatrix.clear();
    computeIndexMaps();
  }
}  

int IndexVector::mIndex(unsigned int n, int m) const {
  return (m+n)/2;
}

void IndexVector::computeIndexMaps() {
  unsigned int i = 0;
  for (int n=0; n<= nmax; n++) {
    for (int n1i=0; n1i <=n; n1i++) {
      cartesian[i] = pair<unsigned int, unsigned int>(n1i,n-n1i);
      cMatrix(n1i, n-n1i) = i;
      polar[i] = pair<unsigned int, int>(n,2*n1i - n);
      pMatrix(n,mIndex(n,2*n1i - n)) = i;
      i++;
    }
  }
}

int IndexVector::getNMax() const {
  return nmax;
}

// size of vector derived from triangular matrix
int IndexVector::getNCoeffs() const {
  return (nmax+1)*(nmax+2)/2;
}

int IndexVector::getN1(unsigned int i) const {
  //return nVector(index,0);
  std::map<unsigned int, std::pair<unsigned int, unsigned int> >::const_iterator cIter = cartesian.find(i);
  if (cIter != cartesian.end())
    return cIter->second.first;
}

int IndexVector::getN2(unsigned int i) const {
  std::map<unsigned int, std::pair<unsigned int, unsigned int> >::const_iterator cIter = cartesian.find(i);
  if (cIter != cartesian.end())
    return cIter->second.second;
}

int IndexVector::getN(unsigned int i) const {
  std::map<unsigned int, std::pair<unsigned int, int> >::const_iterator pIter = polar.find(i);
  if (pIter != polar.end())
    return pIter->second.first;
}

int IndexVector::getM(unsigned int i) const {
  std::map<unsigned int, std::pair<unsigned int, int> >::const_iterator pIter = polar.find(i);
  if (pIter != polar.end())
    return pIter->second.second;
}

unsigned int IndexVector::getCartesianIndex(unsigned int n1, unsigned int n2) const {
  return cMatrix(n1,n2);
}

unsigned int IndexVector::getPolarIndex(unsigned int n, int m) const {
  return pMatrix(n,mIndex(n,m));
}

