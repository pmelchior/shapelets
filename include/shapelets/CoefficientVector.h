#ifndef COEFFICIENTVECTOR_H
#define COEFFICIENTVECTOR_H

#include <NumMatrix.h>
#include <NumVector.h>
#include <Typedef.h>
#include <shapelets/IndexVector.h>

/// Class for storing shapelet coefficients as vector.
/// To perform linear transformation in shapelet space, it is necessary to
///  transform the coefficient matrix into a vector. 
/// In principle, it is completely arbitrary how this should be done, but 
/// in practisesome storage scheme are more convenient than others.
///
/// The vector - matrix mapping is obtained from IndexVector.
///
/// Example:
/// \code
/// // transform coeffs into vector form
///  CoefficientVector<data_t> f(cartesianCoeffs);
///  // perform the convolution
///  CoefficientVector<data_t> h = P*f;
///  // transform vector back to coefficient matrix
///  h.fillCoeffMatrix(cartesianCoeffs);
/// \endcode
/// 
template <class T>
class CoefficientVector : public NumVector<T> {
 public:
  /// Default constructor.
  CoefficientVector(): NumVector<T>() {}

  /// Constructor for CoefficientVector of given \f$n_{max}\f$.
  CoefficientVector(unsigned int nmax): NumVector<T> () {
    CoefficientVector<T>::setNMax(nmax);
  }

  /// Argumented constructor from a NumMatrix.
  /// If the data types of <tt>coeffMatrix</tt> and <tt>*this</tt> differ, the input matrix decides
  /// on the mapping (<tt>data_t</tt> for cartesian, <tt>complex<data_t></tt> for polar coeffs), and
  /// a explicit typecast is performed.
  template <class R> 
    CoefficientVector(const NumMatrix<R>& coeffMatrix) : NumVector<T>() {
    CoefficientVector<T>::setCoeffs(coeffMatrix);
  }

  /// Copy constructor.
  template <class R>
    CoefficientVector(const NumVector<R>& v) {
    if (nVector.getNMax() != computeNMax(v.size()))
      nVector.setNMax(computeNMax(v.size()));
    NumVector<T>::operator=(v);
  } 

  /// Copy constructor.
  template <class R>
    CoefficientVector(const CoefficientVector<R>& cv) {
    if (nVector.getNMax() != cv.nVector.getNMax())
      nVector = cv.nVector;
    NumVector<T>::operator=(cv.getVector());
  }

  const T& operator()(unsigned int i) const {
    return NumVector<T>::operator()(i);
  }
  T& operator()(unsigned int i) {
    return NumVector<T>::operator()(i);
  }

  /// const Matrix-like access operator.
  const T& operator()(int i, int j) const {
    T ownval;
    data_t data_tval;
    complex<data_t> complexval;
    // cartesian coeffs
    if (typeid(ownval) == typeid(data_tval))
      return NumVector<T>::operator()(nVector.getCartesianIndex(i,j));
    // polar coeffs
    else if (typeid(ownval) == typeid(complexval))
      return NumVector<T>::operator()(nVector.getPolarIndex(i,j));
  }

  /// Matrix-like access operator.
  T& operator()(int i, int j) {
    T ownval;
    data_t data_tval;
    complex<data_t> complexval;
    // cartesian coeffs
    if (typeid(ownval) == typeid(data_tval))
      return NumVector<T>::operator()(nVector.getCartesianIndex(i,j));
    // polar coeffs
    else if (typeid(ownval) == typeid(complexval))
      return NumVector<T>::operator()(nVector.getPolarIndex(i,j));
  }

  /// Set new coefficients in matrix form.
  /// If \f$n_{max}\f$ of new coefficient matrix is different, 
  /// CoefficientVector<T> will change accordingly.
  template <class R>
  void setCoeffs(const NumMatrix<R>& coeffMatrix) {
    // set up IndexVector
    int nmax = coeffMatrix.getRows() - 1;
    if (nVector.getNMax() != nmax) {
      nVector.setNMax(nmax);
      // resize storage contained
      NumVector<T>::resize(nVector.getNCoeffs());
    }
    // map matrix onto vector
    // since polar and cartesian coeffs are stored in a different fashion, check datatype first
    R testval;
    T ownval;
    data_t data_tval;
    complex<data_t> complexval;
    // cartesian coeffs
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

  /// Reverse mapping onto a matrix.
  /// If the data types of <tt>coeffMatrix</tt> and <tt>*this</tt> differ, the input matrix decides
  /// on the mapping (<tt>data_t</tt> for cartesian, <tt>complex<data_t></tt> for polar coeffs), and
  /// a explicit typecast is performed.
  template <class R> 
    void fillCoeffMatrix(NumMatrix<R>& coeffMatrix) const {
    // resize matrix if necessary
    int nmax = nVector.getNMax();
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

  /// Get NumMatrix from CoefficientVector.
  /// Mapping is identical to fillCoeffMatrix().
  NumMatrix<T> getCoeffMatrix() const {
    NumMatrix<T> M;
    CoefficientVector<T>::fillCoeffMatrix(M);
    return M;
  }

  /// Set \f$n_{max}\f$.
  void setNMax(unsigned int nmax) {
    if (nVector.getNMax() != nmax) {
      nVector.setNMax(nmax);
      NumVector<T>::resize(nVector.getNCoeffs());
      NumVector<T>::clear();
    }
  }

  /// Get \f$n_{max}\f$.
  unsigned int getNMax() const {
    return nVector.getNMax();
  }
  /// Get \f$n_{coeffs}\f$.
  unsigned int getNCoeffs() const {
    return NumVector<T>::size();
  }

  /// Get the IndexVector vector <-> matrix mapping.
  const IndexVector& getIndexVector() const {
    return nVector;
  }
  const NumVector<T>& getNumVector() const {
    return *this;
  }
 private:
  IndexVector nVector;
  int computeNMax(int nCoeffs) {
    return (int)floor(-1.5 + 0.5*sqrt(9.0 + 8*(nCoeffs-1))); 
  }
};

template <class T>
NumVector<T> operator*(const NumMatrix<T>& M, const CoefficientVector<T>& cv) {
  return M*cv.getNumVector();
}
#endif
