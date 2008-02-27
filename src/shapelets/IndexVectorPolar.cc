#include <shapelets/IndexVectorPolar.h>

using namespace std;

IndexVectorPolar::IndexVectorPolar () {
  nmax = 0;
}

IndexVectorPolar::IndexVectorPolar (unsigned int nmax) {
  setNMax(nmax);
}

void IndexVectorPolar::setNMax(int innmax) {
  if (innmax != nmax) {
    nmax = innmax;
    polar.clear();
    pMatrix.resize(nmax+1,nmax+1);
    pMatrix.clear();
    computeIndexMaps();
  }
}  

unsigned int IndexVectorPolar::getNMax() const {
  return nmax;
}

unsigned int IndexVectorPolar::getNCoeffs() const {
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
}

int IndexVectorPolar::getState2(unsigned int i) const {
  std::map<unsigned int, std::pair<unsigned int, int> >::const_iterator pIter = polar.find(i);
  if (pIter != polar.end())
    return pIter->second.second;
}

unsigned int IndexVectorPolar::getIndex1(unsigned int i) const {
  return IndexVectorPolar::getState1(i);
}

unsigned int IndexVectorPolar::getIndex2(unsigned int i) const {
  std::map<unsigned int, std::pair<unsigned int, int> >::const_iterator pIter = polar.find(i);
  if (pIter != polar.end())
    return mIndex(pIter->second.first,pIter->second.second);
}

unsigned int IndexVectorPolar::getIndex(int n, int m) const {
  return pMatrix(n,mIndex(n,m));
}

