#include <MatrixManipulations.h>
#include <math.h>
#include <iostream>

namespace ublas = boost::numeric::ublas;
typedef complex<double> Complex;

// maps the range m = -n, -n +2, .. , n to
// m = 0,1,..,n to efficiently store in matrix
unsigned int mIndex(int m, int n) {
  return (m + n)/2;
}

// mapping a triangular matrix onto a vector of minimal size
// since the matrices have either lower left (polar) or upper left form
// the counting is done in different fashion
void matrixMapping(const NumMatrix<double>& matrix, NumVector<double>& vector,bool polar, const IndexVector& nVector) {
  if (vector.size() != nVector.getNCoeffs())
    vector.resize(nVector.getNCoeffs());
  for (int n =0; n< nVector.getNCoeffs(); n++ ){
    if (polar) vector(n) = matrix(nVector.getN(n),nVector.getN1(n));
    else vector(n) = matrix(nVector.getN1(n),nVector.getN2(n));
  }
}

void matrixMapping(const NumMatrix<Complex>& matrix, NumVector<Complex>& vector,bool polar, const IndexVector& nVector) {
  if (vector.size() != nVector.getNCoeffs())
    vector.resize(nVector.getNCoeffs());
  for (int n =0; n< nVector.getNCoeffs(); n++ ){
    if (polar) vector(n) = matrix(nVector.getN(n),nVector.getN1(n));
    else vector(n) = matrix(nVector.getN1(n),nVector.getN2(n));
  }
}

void matrixMapping(const NumMatrix<double>& matrix, NumVector<Complex>& vector,bool polar, const IndexVector& nVector) {
  if (vector.size() != nVector.getNCoeffs())
    vector.resize(nVector.getNCoeffs());
  for (int n =0; n< nVector.getNCoeffs(); n++ ){
    if (polar) vector(n) = matrix(nVector.getN(n),nVector.getN1(n));
    else vector(n) = matrix(nVector.getN1(n),nVector.getN2(n));
  }
}

void matrixMapping(const NumMatrix<double>& matrix,ublas::matrix_row< ublas::matrix<double> >& mr,bool polar, const IndexVector& nVector) {
  for (int n =0; n< nVector.getNCoeffs(); n++ ){
    if (polar) mr(n) = matrix(nVector.getN(n),nVector.getN1(n));
    else mr(n) = matrix(nVector.getN1(n),nVector.getN2(n));
  }
}

// the reverse operation to matrixMapping
// here for the Complex values of the polar vector/matrix
void vectorMapping(const ublas::vector<Complex>& vector,NumMatrix<Complex>& matrix, const IndexVector& nVector) {
  for (int n =0; n< nVector.getNCoeffs(); n++ ){
    matrix(nVector.getN(n),nVector.getN1(n)) = vector(n);
  }
}

// the reverse operation to matrixMapping
// for double vector -> double matrix
void vectorMapping(const ublas::vector<double>& vector,NumMatrix<double>& matrix, const IndexVector& nVector) {
  for (int n =0; n< nVector.getNCoeffs(); n++ ){
    matrix(nVector.getN1(n),nVector.getN2(n)) = vector(n);
  }
}

// the reverse operation to matrixMapping
// here for the double values of the cartesian vector/matrix
// when coming from a polar transformation the matrix is Complex
void vectorMapping(const ublas::vector<Complex>& vector,NumMatrix<double>& matrix, const IndexVector& nVector) {
  for (int n =0; n< nVector.getNCoeffs(); n++ ){
    // this is a empirical value for the transformations done with coeffs
    if (fabs(imag(vector(n))) > 0.01*fabs(real(vector(n))) && fabs(imag(vector(n))) > 1e-10) 
      std::cout << "# Warning: Polar -> Cartesian Transformation: polar matrix describes complex function! " << n << " " << vector(n) << std::endl;
    matrix(nVector.getN1(n),nVector.getN2(n)) = real(vector(n));
  }
}
