#include <MatrixManipulations.h>
#include <math.h>
#include <iostream>

namespace ublas = boost::numeric::ublas;
typedef complex<double> Complex;

// size of vector derived from triangular matrix
int getNCoeffs(int nmax) {
  return (nmax+1)*(nmax+2)/2;
}

int getNMax(int nCoeffs) {
  return (int)floor(-1.5 + 0.5*sqrt(9.0 + 8*(nCoeffs-1)));
}
// maps the range m = -n, -n +2, .. , n to
// m = 0,1,..,n to efficiently store in matrix
unsigned int mIndex(int m, int n) {
  return (m + n)/2;
}
// mapping of the vector index to the matrix indices
int getN1(const NumMatrix<int>& nVector,int i) {
  return nVector(0,i);
}
int getN2 (const NumMatrix<int>& nVector,int i) {
  return nVector(1,i);
}
int getN1Polar(const NumMatrix<int>& nVector,int i) {
  return nVector(2,i);
}

// build coefficient mapping vector
// this is voodoo!
// it creates vectors which map the diagonal elements
// of a matrix in a row for accessing them in a vector
void makeNVector(NumMatrix<int>& nVector,int nCoeffs, int nmax) {
  if (nVector.getRows() != 3 || nVector.getColumns() != nCoeffs)
    nVector = NumMatrix<int>(3,nCoeffs);
  int i=0;
  for (int n=0; n<= nmax; n++) {
    for (int n1i=0; n1i <=n; n1i++) {
      // the first two map diagonal elements
      nVector(0,i) = n1i;
      nVector(1,i) = n-n1i;
      // this one is needed for the polar matrix representation
      nVector(2,i) = n;
      i++;
    }
  }
}

// mapping a triangular matrix onto a vector of minimal size
// since the matrices have either lower left (polar) or upper left form
// the counting is done in different fashion
void matrixMapping(const NumMatrix<double>& matrix, NumVector<double>& vector,bool polar, const NumMatrix<int>& nVector,int nCoeffs) {
//void matrixMapping(NumMatrix<double>& matrix,NumVector<double>& vector,bool polar,NumMatrix<int>& nVector,int nCoeffs)
  for (int n =0; n< nCoeffs; n++ ){
    if (polar) vector(n) = matrix(getN1Polar(nVector,n),getN1(nVector,n));
    else vector(n) = matrix(getN1(nVector,n),getN2(nVector,n));
  }
}

void matrixMapping(const NumMatrix<Complex>& matrix, NumVector<Complex>& vector,bool polar, const NumMatrix<int>& nVector, int nCoeffs) {
  for (int n =0; n< nCoeffs; n++ ){
    if (polar) vector(n) = matrix(getN1Polar(nVector,n),getN1(nVector,n));
    else vector(n) = matrix(getN1(nVector,n),getN2(nVector,n));
  }
}

void matrixMapping(const NumMatrix<double>& matrix, NumVector<Complex>& vector,bool polar, const NumMatrix<int>& nVector, int nCoeffs) {
  for (int n =0; n< nCoeffs; n++ ){
    if (polar) vector(n) = matrix(getN1Polar(nVector,n),getN1(nVector,n));
    else vector(n) = matrix(getN1(nVector,n),getN2(nVector,n));
  }
}

void matrixMapping(const NumMatrix<double>& matrix,ublas::matrix_row< ublas::matrix<double> >& mr,bool polar, const NumMatrix<int>& nVector,int nCoeffs) {
  for (int n =0; n< nCoeffs; n++ ){
    if (polar) mr(n) = matrix(getN1Polar(nVector,n),getN1(nVector,n));
    else mr(n) = matrix(getN1(nVector,n),getN2(nVector,n));
  }
}

// the reverse operation to matrixMapping
// here for the Complex values of the polar vector/matrix
void vectorMapping(const ublas::vector<Complex>& vector,NumMatrix<Complex>& matrix, const NumMatrix<int>& nVector,int nCoeffs) {
  for (int n =0; n< nCoeffs; n++ ){
    matrix(getN1Polar(nVector,n),getN1(nVector,n)) = vector(n);
  }
}

// the reverse operation to matrixMapping
// for double vector -> double matrix
void vectorMapping(const ublas::vector<double>& vector,NumMatrix<double>& matrix, const NumMatrix<int>& nVector,int nCoeffs) {
  for (int n =0; n< nCoeffs; n++ ){
    matrix(getN1(nVector,n),getN2(nVector,n)) = vector(n);
  }
}

// the reverse operation to matrixMapping
// here for the double values of the cartesian vector/matrix
// when coming from a polar transformation the matrix is Complex
void vectorMapping(const ublas::vector<Complex>& vector,NumMatrix<double>& matrix, const NumMatrix<int>& nVector,int nCoeffs) {
  for (int n =0; n< nCoeffs; n++ ){
    // this is a empirical value for the transformations done with coeffs
    if (fabs(imag(vector(n))) > 0.01*fabs(real(vector(n))) && fabs(imag(vector(n))) > 1e-10) 
      std::cout << "# Warning: Polar -> Cartesian Transformation: polar matrix describes complex function! " << n << " " << vector(n) << std::endl;
    matrix(getN1(nVector,n),getN2(nVector,n)) = real(vector(n));
  }
}
