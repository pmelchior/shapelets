#ifndef COEFFICIENTVECTOR_H
#define COEFFICIENTVECTOR_H

#include <NumMatrix.h>
#include <NumVector.h>
#include <Typedef.h>
#include <shapelets/IndexVector.h>

/// Class for storing shapelet coefficients as vector.
/// To perform linear transformation in shapelet sace, it is necessary to transform the coefficient 
/// matrix into a vector. In principle, it is completely arbitrary how this should be done, but in practise
/// some stroage scheme are more convenient than others.
///
/// In this shapelet implementation the matrix storage scheme for cartesian and polar coefficient differs, 
/// and so the mapping onto a coefficient vector (accomplished by IndexVector) also has to differ:
/// - for cartesian coefficients:
/// \f[
/// \begin{bmatrix}
/// (0,0) & (0,1) & (0,2) & (0,3)\\
/// (1,0) & (1,1) & (1,2) & \\
/// (2,0) & (2,1) & & \\
/// (3,0) & & &
/// \end{bmatrix}
/// \longmapsto
/// \bigl[ (0,0), (0,1), (1,0), (0,2), (1,1), (2,0), (0,3), (1,2), (2,1), (3,0)\bigr]
/// \f]
/// - for polar coefficients:
/// \f[
/// \begin{bmatrix}
/// (0,0) & & & \\
/// (1,-1) & (1,1) & & \\
/// (2,-2) & (2,0) & (2,2) & \\
/// (3,-3) & (3,-1) & (3,1) & (3,3)
/// \end{bmatrix}
/// \longmapsto
/// \bigl[ (0,0), (1,-1), (1,1), (2,-2), (2,0), (2,2), (3,-3), (3,-1), (3,1), (3,3)\bigr]
/// \f]
///
/// Usage is straightforward:
/// \code
/// // transform coeffs into vector form
///  CoefficientVector<data_t> f(cartesianCoeffs);
///  // perform the convolution
///  CoefficientVector<data_t> h = P*f;
///  // transform vector back to coefficient matrix
///  h.fillCoefficientMatrix(cartesianCoeffs);
/// \endcode
/// 
template <class T>
class CoefficientVector : public NumVector<T> {
 public:
  /// Default constructor.
  CoefficientVector(): NumVector<T>() {
  }
  CoefficientVector(unsigned int size): NumVector<T> (size) {
    nVector = IndexVector(computeNMax(size));
  }
  /// Argumented constructor from a NumMatrix.
  /// If the data types of <tt>coeffMatrix</tt> and <tt>*this</tt> differ, the input matrix decides
  /// on the mapping (<tt>data_t</tt> for cartesian, <tt>complex<data_t></tt> for polar coeffs), and
  /// a explicit typecast is performed.
  template <class R> 
    CoefficientVector(const NumMatrix<R>& coeffMatrix) : NumVector<T>() {
    // set up IndexVector
    int nmax = coeffMatrix.getRows() - 1;
    nVector = IndexVector(nmax);
    // resize storage contained
    NumVector<T>::resize(nVector.getNCoeffs());
    // map matrix onto vector
    // since polar and cartesian coeffs are stored in a different fashion, check datatype first
    R testval;
    T ownval;
    data_t data_tval;
    complex<data_t> complexval;
    // data_t cartesian coeffs
    if (typeid(testval) == typeid(data_tval)) {
      for (int n = 0; n < NumVector<T>::size(); n++) {
	// no typecast required
	if (typeid(testval) == typeid(ownval))
	  NumVector<T>::operator()(n) = coeffMatrix(nVector.getN1(n),nVector.getN2(n));
	// data_t -> complex typecast
	else
	  NumVector<T>::operator()(n) = static_cast<T>(coeffMatrix(nVector.getN1(n),nVector.getN2(n)));
      }
    } 
    // polar coeffs
    else if (typeid(testval) == typeid(complexval)) {
      for (int n = 0; n < NumVector<T>::size(); n++) {
	// no typecast required
	if (typeid(testval) == typeid(ownval))
	  NumVector<T>::operator()(n) = coeffMatrix(nVector.getN(n),nVector.getN1(n));
	// complex -> data_t typecast
	else
	  NumVector<T>::operator()(n) = static_cast<T>(coeffMatrix(nVector.getN1(n),nVector.getN2(n)));
      }
    }
  }
  /// Copy constructor.
  template <class R>
    CoefficientVector(const NumVector<R>& v) {
    NumVector<T>::operator=(v);
    nVector = IndexVector(computeNMax(v.size()));
  } 
  /// Copy constructor.
  template <class R>
    CoefficientVector(const CoefficientVector<T>& cv) {
    NumVector<T>::operator=(cv);
    nVector = cv.nVector;
  }
  /// Reverse mapping onto a matrix.
  /// If the data types of <tt>coeffMatrix</tt> and <tt>*this</tt> differ, the input matrix decides
  /// on the mapping (<tt>data_t</tt> for cartesian, <tt>complex<data_t></tt> for polar coeffs), and
  /// a explicit typecast is performed.
  template <class R> 
    void fillCoeffMatrix(NumMatrix<R>& coeffMatrix) const {
    // resize matrix if necessary
    int nmax = nVector.getOrder();
    if (coeffMatrix.getRows() != nmax + 1 || coeffMatrix.getColumns() != nmax +1)
      coeffMatrix.resize_clear(nmax+1,nmax+1);
    // map vector onto matrix
    // since polar and cartesian coeffs are stored in a different fashion, check datatype first
    R testval;
    T ownval;
    data_t data_tval;
    complex<data_t> complexval;
    // cartesian coeffs
    if (typeid(testval) == typeid(data_tval)) {
      for (int n = 0; n < NumVector<T>::size(); n++) {
	// no typecast
	if (typeid(testval) == typeid(ownval))
	  coeffMatrix(nVector.getN1(n),nVector.getN2(n)) = NumVector<T>::operator()(n);
	// data_t -> complex typecast
	else
	  coeffMatrix(nVector.getN1(n),nVector.getN2(n)) = static_cast<T>(NumVector<T>::operator()(n));
      }
    } 
    // polar coeffs
    else if (typeid(testval) == typeid(complexval)) {
      for (int n = 0; n < NumVector<T>::size(); n++) {
	// no typecast
	if (typeid(testval) == typeid(ownval))
	  coeffMatrix(nVector.getN(n),nVector.getN1(n)) = NumVector<T>::operator()(n);
	// complex -> data_t typecast
	else
	  coeffMatrix(nVector.getN(n),nVector.getN1(n)) = static_cast<T>(NumVector<T>::operator()(n));
      }
    }
  }
  /// Get the IndexVector vector <-> matrix mapping.
  const IndexVector& getIndexVector() const {
    return nVector;
  }
 private:
  IndexVector nVector;
  int computeNMax(int nCoeffs) {
    return (int)floor(-1.5 + 0.5*sqrt(9.0 + 8*(nCoeffs-1))); 
  }
};

#endif
