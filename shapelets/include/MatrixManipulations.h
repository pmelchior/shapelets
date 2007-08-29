#ifndef MATRIXMANIPULATIONS_H
#define MATRIXMANIPULATIONS_H

#include <NumMatrix.h>
#include <NumVector.h>
#include <IndexVector.h>
#include <complex.h>
#include <iostream>

/// Functions for manipulating matrices.
/// When working matrices of coefficients, it's often more convenient
/// to transform the tensor into a matrix and the matrix into a vector
/// of appropriate size.

/// Maps the range m = -n, -n +2, .. , n to
/// m = 0,1,..,n to efficiently store in matrix.
unsigned int mIndex(int m, int n);

/// Map a triangular matrix onto a vector of minimal size.
/// Since the matrices have either lower left (polar) or upper left form,
/// the counting is done in different fashion
void matrixMapping(const NumMatrix<double>& matrix, NumVector<double>& vector,bool polar, const IndexVector& nVector);
/// Map a triangular matrix onto a vector of minimal size.
/// Counting in polar or cartesian fashion.
void matrixMapping(const NumMatrix< complex<double> >& matrix, NumVector< complex<double> >& vector,bool polar, const IndexVector& nVector);
/// Map a triangular matrix onto a vector of minimal size.
/// Counting in polar or cartesian fashion.
void matrixMapping(const NumMatrix<double>& matrix, NumVector< complex<double> >& vector,bool polar, const IndexVector& nVector);
/// Map a triangular matrix onto a matrix_row of minimal size.
/// Counting in polar or cartesian fashion.
void matrixMapping(const NumMatrix<double>& matrix,boost::numeric::ublas::matrix_row< boost::numeric::ublas::matrix<double> >& mr,bool polar, const IndexVector& nVector);

/// Reverse operation to matrixMapping().
/// Here for the Complex values of the polar vector/matrix.
void vectorMapping(const boost::numeric::ublas::vector< complex<double> >& vector,NumMatrix< complex<double> >& matrix, const IndexVector& nVector);
/// Reverse operation to matrixMapping().
/// Here for double vector -> double matrix.
void vectorMapping(const boost::numeric::ublas::vector<double>& vector,NumMatrix<double>& matrix, const IndexVector& nVector);
/// Reverse operation to matrixMapping().
/// Here for the double values of the cartesian vector/matrix 
/// (when coming from a polar transformation the matrix is Complex)
void vectorMapping(const boost::numeric::ublas::vector< complex<double> >& vector,NumMatrix<double>& matrix, const IndexVector& nVector);


/// Transform arbitrary matrix into triangular matrix of appropraite dimension.
/// Copies entries from input matrix into lower left corner of the triangular matrix
/// as long as there are entries unequal to 0 on the diagonal. 
template <class T>
int triangularizeCoeffs(const NumMatrix<T>& input, NumMatrix<T>& triangular, double cutoff) {
  int nmax = getNMax(input,cutoff);
  triangular = NumMatrix<T>(nmax+1,nmax+1);
  for (int i = 0; i< nmax +1; i++) 
     for (int j = 0; j< nmax+1; j++) 
       // truncation at nmax
       if (i < input.getRows() && j < input.getColumns() && i+j <= nmax)
	 triangular(i,j) = input(i,j);
  return nmax;
}

/// Get maximal order to keep for the triangular matrix.
/// Reads the matrix from top right to bottom left along the diagonals and returns
/// the number of the diagonal where not all entries are 0.
template <class T>
int getNMax(const NumMatrix<T>& matrix, double cutoff) {
//   // first determine the maximum entry
//   double max = 0;
//   for (int i = 0; i < coeffs.getRows(); i++)
//     for (int j = 0; j < coeffs.getColumns(); j++)
//       if (fabs(coeffs(i,j)) > max) max = fabs(coeffs(i,j));
  
//   // FIXME: should be fine-tuned by decomposing a catalog of galaxies
//   // should be also a globally defined value
//   double cutoff = 1e-5*max;
//   cout << "# Truncation: cutoff = " << cutoff << endl;

  int result;
  // start in lower right corner and work back to top left corner
  for (int diag = (matrix.getRows() + matrix.getColumns() -2); diag>=0; diag-- ) {
    result = diag;
    bool isempty = 1;
    int startj, stopj;
    if (diag >= matrix.getRows()-1) {
      startj = diag - matrix.getRows() + 1;
      stopj = matrix.getColumns() - 1;
    } else {
      startj = 0;
      stopj = diag;
    }
    for (int j = startj; j <= stopj; j++) {
      int i = diag -j;
       if (abs(matrix(i,j)) > cutoff) {
	isempty = 0;
	break;
      }
    }
    if (!isempty) break;
  }
  return result;
}


#endif
