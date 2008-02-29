#ifndef INDEXVECTOR_H
#define INDEXVECTOR_H

#include <NumMatrix.h>
#include <map>

/// Index of vector representation of coefficient matrix.
/// It is often convenient to store the shapelet coefficients (matrix) as vectors.
/// This class provides the methods to infer the matrix indices or the 
/// eigenstate numbers from the vector index and vice versa.\n\n
/// The storage scheme employs \em triangular layout, because this allows
/// direct conversion between cartesian and polar states (cf. PolarTransformation).
/// The \em triangular mapping looks as follows:
/// - for cartesian coefficients:
/// \f[ \begin{bmatrix} (0,0) & (0,1) & (0,2) & (0,3)\\ (1,0) & (1,1) & (1,2) & \\ (2,0) & (2,1) & & \\ (3,0) & & &  \end{bmatrix}  \longmapsto  \bigl[ (0,0), (0,1), (1,0), (0,2), (1,1), (2,0), (0,3), (1,2), (2,1), (3,0)\bigr]  \f]
/// - for polar coefficients:
/// \f[ \begin{bmatrix} (0,0) & & & \\ (1,-1) & (1,1) & & \\ (2,-2) & (2,0) & (2,2) & \\ (3,-3) & (3,-1) & (3,1) & (3,3) \end{bmatrix} \longmapsto \bigl[ (0,0), (1,-1), (1,1), (2,-2), (2,0), (2,2), (3,-3), (3,-1), (3,1), (3,3)\bigr] \f]
///
/// Other storage schemes can be implemented here easily.

class IndexVector {
 public:
  virtual ~IndexVector() {};
  /// Set \f$n_{max}\f$.
  virtual void setNMax(unsigned int nmax) = 0;
  /// Get \f$n_{max}\f$.
  virtual unsigned int getNMax() const = 0;
  /// Get \f$n_{coeffs}\f$.
  virtual unsigned int getNCoeffs() const = 0;
  /// Get the first eigenstate number from the vector index.
  virtual int getState1(unsigned int index) const = 0;
  /// Get the second eigenstate number from the vector index.
  virtual int getState2(unsigned int index) const = 0;
  /// Get the first matrix index from vector index.
  virtual unsigned int getIndex1(unsigned int index) const = 0;
  /// Get the second matrix index from vector index.
  virtual unsigned int getIndex2(unsigned int index) const = 0;
  /// Get the vector index from the eigenstate numbers.
  virtual unsigned int getIndex(int i, int j) const = 0;
};

#endif
