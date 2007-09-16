#ifndef INDEXVECTOR_H
#define INDEXVECTOR_H

#include <NumMatrix.h>

/// Index of vector representation of coefficient matrix.
/// It is often convenient to store the shapelet coefficients (matrix) as a vector.
/// This class provides the methods to recover the matrix indices from the vector index.

class IndexVector : public NumMatrix<int> {
 public:
  /// Default constructor.
  IndexVector();
  /// Argumented constructor.
  IndexVector (int nmax);
  /// Set \f$n_{max}\f$ of the associated coefficient matrix.
  void setOrder(int nmax);
  /// Set \f$n_{max}\f$ of the associated coefficient matrix.
  int getOrder() const;
  /// Get the number of shapelet coefficients.
  int getNCoeffs() const;
  /// Get the index \f$n_1\f$ of the associated cartesian coefficient matrix from the vector index.
  int getN1(unsigned int index) const;
  /// Get the index \f$n_2\f$ of the associated cartesian coefficient matrix from the vector index.
  int getN2(unsigned int index) const;
  /// Get the index \f$n\f$ of the associated polar coefficient matrix from the vector index.
  int getN(unsigned int index) const;
  /// Get the index \f$m\f$ of the associated polar coefficient matrix from the vector index.
  int getM(unsigned int index) const;
 private:
  void computeIndexVector();
  unsigned int nmax;
  NumMatrix<int> nVector;
};
#endif
