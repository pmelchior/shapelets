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
  /// Default constructor.
  IndexVector();
  /// Argumented constructor.
  IndexVector (int nmax);
  /// Set \f$n_{max}\f$ of the associated coefficient matrix.
  void setNMax(int nmax);
  /// Set \f$n_{max}\f$ of the associated coefficient matrix.
  int getNMax() const;
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
  /// Get the vector index from the cartesian eigenstates numbers.
  unsigned int getCartesianIndex(unsigned int n1, unsigned int n2) const;
  /// Get the vector index from the polar eigenstate numbers.
  unsigned int getPolarIndex(unsigned int n, int m) const;
 private:
  void computeIndexMaps();
  unsigned int nmax;
  std::map<unsigned int, std::pair<unsigned int, unsigned int> > cartesian;
  std::map<unsigned int, std::pair<unsigned int, int> > polar;
  NumMatrix<unsigned int> cMatrix;
  NumMatrix<unsigned int> pMatrix;
  int mIndex(unsigned int n, int m) const;
};
#endif
