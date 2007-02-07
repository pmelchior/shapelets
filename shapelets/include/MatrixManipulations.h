#ifndef MATRIXMANIPULATIONS_H
#define MATRIXMANIPULATIONS_H

#include <NumMatrix.h>
#include <NumVector.h>
//#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
//#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <complex.h>
#include <iostream>

/// Functions for manipulating matrices.
/// When working matrices of coefficients, it's often more convenient
/// to transform the tensor into a matrix and the matrix into a vector
/// of appropriate size.
/// Unfortunately templated functions don't work here, so we have to set
/// the type for the matrices explicitly.
/// \todo nVector should be own class.

/// Size of vector derived from triangular matrix.
int getNCoeffs(int nmax);

/// Maximum shapelets order from number of coefficients
int getNMax(int nCoeffs);

/// Maps the range m = -n, -n +2, .. , n to
/// m = 0,1,..,n to efficiently store in matrix.
unsigned int mIndex(int m, int n);

/// Mapping of vector index to matrix indices.
/// The line of a cartesian matrix.
int getN1(const NumMatrix<int>& nVector,int i);

/// Mapping of vector index to matrix indices.
/// The column of the matrix.
int getN2 (const NumMatrix<int>& nVector,int i);

/// Mapping of vector index to matrix indices.
/// The line of a polar matrix.
int getN1Polar(const NumMatrix<int>& nVector,int i);

/// Build polar coefficient mapping vector.
/// It creates vectors which map the diagonal elements
/// of a matrix in one row for accessing them successively in a vector
void makeNVector(NumMatrix<int>& nVector,int nCoeffs, int nmax);

/// Map a triangular matrix onto a vector of minimal size.
/// Since the matrices have either lower left (polar) or upper left form,
/// the counting is done in different fashion
void matrixMapping(const NumMatrix<double>& matrix, NumVector<double>& vector,bool polar, const NumMatrix<int>& nVector,int nCoeffs);
/// Map a triangular matrix onto a vector of minimal size.
/// Counting in polar or cartesian fashion.
void matrixMapping(const NumMatrix< complex<double> >& matrix, NumVector< complex<double> >& vector,bool polar, const NumMatrix<int>& nVector,int nCoeffs);
/// Map a triangular matrix onto a vector of minimal size.
/// Counting in polar or cartesian fashion.
void matrixMapping(const NumMatrix<double>& matrix, NumVector< complex<double> >& vector,bool polar, const NumMatrix<int>& nVector,int nCoeffs);
/// Map a triangular matrix onto a matrix_row of minimal size.
/// Counting in polar or cartesian fashion.
void matrixMapping(const NumMatrix<double>& matrix,boost::numeric::ublas::matrix_row< boost::numeric::ublas::matrix<double> >& mr,bool polar, const NumMatrix<int>& nVector,int nCoeffs);

/// Reverse operation to matrixMapping().
/// Here for the Complex values of the polar vector/matrix.
void vectorMapping(const boost::numeric::ublas::vector< complex<double> >& vector,NumMatrix< complex<double> >& matrix, const NumMatrix<int>& nVector,int nCoeffs);
/// Reverse operation to matrixMapping().
/// Here for double vector -> double matrix.
void vectorMapping(const boost::numeric::ublas::vector<double>& vector,NumMatrix<double>& matrix, const NumMatrix<int>& nVector,int nCoeffs);
/// Reverse operation to matrixMapping().
/// Here for the double values of the cartesian vector/matrix 
/// (when coming from a polar transformation the matrix is Complex)
void vectorMapping(const boost::numeric::ublas::vector< complex<double> >& vector,NumMatrix<double>& matrix, const NumMatrix<int>& nVector,int nCoeffs);


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

/// Change the dimension of a square matrix by \f$\delta\f$ in every dimension.
/// \f$\delta\f$ can be positive (adding zeros) or negative (slicing).
template <class T>
void changeDimension(NumMatrix<T>& matrix, int delta) {
  int oldDim = matrix.getRows();
  boost::numeric::ublas::matrix<T> tmpMatrix(oldDim + delta, oldDim + delta);
  if (delta > 0) {
    boost::numeric::ublas::matrix_range<boost::numeric::ublas::matrix<T> > (tmpMatrix, boost::numeric::ublas::range (0, oldDim), boost::numeric::ublas::range (0, oldDim)) = matrix.getMatrix();
    matrix.getMatrix() = tmpMatrix;
  } else if (delta < 0) {
    tmpMatrix = boost::numeric::ublas::matrix_range<boost::numeric::ublas::matrix<T> > (matrix.getMatrix(), boost::numeric::ublas::range (0, oldDim + delta), boost::numeric::ublas::range (0, oldDim + delta));
    matrix.getMatrix() = tmpMatrix;
  }
}

/// Transform vector of \f$ n\cdot m\f$ entries into \f$ n\times m\f$ matrix.
/// If columnwise == 1, the matrix will be a \f$ m\times n\f$ one.
template <class T>
void transformVector2Matrix(NumVector<T>& V, NumMatrix<T>& M, unsigned int n, bool columnwise) {
  if (V.size()%n==0) {
    int m = V.size()/n;
    for (int i=0; i<m; i++){
      if (columnwise) {
	if (M.getRows() != m || M.getColumns() != n)
	  M = NumMatrix<T>(m,n);
	boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<T> > (M.getMatrix(),i) =
	  boost::numeric::ublas::vector_range<boost::numeric::ublas::vector<T> > (V, boost::numeric::ublas::range(i*n,(i+1)*n));
      } else {
	if (M.getRows() != n || M.getColumns() != m)
	  M = NumMatrix<T>(n,m);
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<T> > (M.getMatrix(),i) =
	  boost::numeric::ublas::vector_range<boost::numeric::ublas::vector<T> > (V, boost::numeric::ublas::range(i*n,(i+1)*n));
      }
    }
  }
}

#endif
