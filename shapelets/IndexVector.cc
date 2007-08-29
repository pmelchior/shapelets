#include <IndexVector.h>

IndexVector::IndexVector () : NumMatrix<int>() {
}

IndexVector::IndexVector (int nmax) : NumMatrix<int>() {
  setOrder(nmax);
}

void IndexVector::setOrder(int innmax) {
  nmax = innmax;
  // set size to 3 x nCoeffs
  if (NumMatrix<int>::getRows() != 3 && NumMatrix<int>::getColumns() != getNCoeffs()) 
    NumMatrix<int>::NumMatrix<int>(3,getNCoeffs());
  computeIndexVector();
}  

void IndexVector::computeIndexVector() {
  int i=0;
  for (int n=0; n<= nmax; n++) {
    for (int n1i=0; n1i <=n; n1i++) {
      // the first two map diagonal elements
      NumMatrix<int>::operator()(0,i) = n1i;
      NumMatrix<int>::operator()(1,i) = n-n1i;
      // this one is needed for the polar matrix representation
      NumMatrix<int>::operator()(2,i) = n;
      i++;
    }
  }
}

int IndexVector::getOrder() const {
  return nmax;
}

// size of vector derived from triangular matrix
int IndexVector::getNCoeffs() const {
  return (nmax+1)*(nmax+2)/2;
}

int IndexVector::getN1(unsigned int i) const {
  return NumMatrix<int>::operator()(0,i);
}

int IndexVector::getN2(unsigned int i) const {
  return NumMatrix<int>::operator()(1,i);
}

int IndexVector::getN(unsigned int i) const {
  return NumMatrix<int>::operator()(2,i);
}
