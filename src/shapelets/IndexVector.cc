#include <shapelets/IndexVector.h>

IndexVector::IndexVector () : NumMatrix<int>() {
}

IndexVector::IndexVector (int nmax) : NumMatrix<int>() {
  setOrder(nmax);
}

void IndexVector::setOrder(int innmax) {
  nmax = innmax;
  // set size to 3 x nCoeffs
  if (NumMatrix<int>::getRows() != getNCoeffs() || NumMatrix<int>::getColumns() != 3) 
    NumMatrix<int>::resize(getNCoeffs(),3);
  computeIndexVector();
}  

void IndexVector::computeIndexVector() {
  int i=0;
  for (int n=0; n<= nmax; n++) {
    for (int n1i=0; n1i <=n; n1i++) {
      // the first two map diagonal elements
      NumMatrix<int>::operator()(i,0) = n1i;
      NumMatrix<int>::operator()(i,1) = n-n1i;
      // this one is needed for the polar matrix representation
      NumMatrix<int>::operator()(i,2) = n;
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

int IndexVector::getN1(unsigned int index) const {
  return NumMatrix<int>::operator()(index,0);
}

int IndexVector::getN2(unsigned int index) const {
  return NumMatrix<int>::operator()(index,1);
}

int IndexVector::getN(unsigned int index) const {
  return NumMatrix<int>::operator()(index,2);
}

int IndexVector::getM(unsigned int index) const {
  int i = getN(index);
  int j = getN1(index);
  return 2*j - i;
}
