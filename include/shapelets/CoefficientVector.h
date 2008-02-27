#ifndef COEFFICIENTVECTOR_H
#define COEFFICIENTVECTOR_H

#include <NumMatrix.h>
#include <NumVector.h>
#include <Typedef.h>
#include <shapelets/IndexVector.h>
#include <shapelets/IndexVectorCartesian.h>
#include <shapelets/IndexVectorPolar.h>

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
  CoefficientVector(): NumVector<T>() {
    if (isComplex()) 
      nVector = new IndexVectorPolar();
    else
      nVector = new IndexVectorCartesian();
  }

  /// Constructor for CoefficientVector of given \f$n_{max}\f$.
  CoefficientVector(unsigned int nmax): NumVector<T> () {
    if (isComplex()) 
      nVector = new IndexVectorPolar();
    else
      nVector = new IndexVectorCartesian();
    CoefficientVector<T>::setNMax(nmax);
  }

  /// Argumented constructor from a NumMatrix.
  CoefficientVector(const NumMatrix<T>& coeffMatrix) : NumVector<T>() {
    if (isComplex()) 
      nVector = new IndexVectorPolar();
    else
      nVector = new IndexVectorCartesian();
    CoefficientVector<T>::setCoeffs(coeffMatrix);
  }

  /// Copy constructor.
  CoefficientVector(const NumVector<T>& v) {
    if (isComplex()) 
      nVector = new IndexVectorPolar(computeNMax(v.size()));
    else
      nVector = new IndexVectorCartesian(computeNMax(v.size()));
    NumVector<T>::operator=(v);
  } 

  /// Copy constructor.
  CoefficientVector(const CoefficientVector<T>& cv) {
    if (isComplex()) 
      nVector = new IndexVectorPolar(cv.nVector->getNMax());
    else
      nVector = new IndexVectorCartesian(cv.nVector->getNMax());
    NumVector<T>::operator=(cv.getNumVector());
  }

  /// Destructor.
  ~CoefficientVector() {
    delete nVector;
  }
  
  /// Copy operator.
  void operator= (const NumVector<T>& v) {
    if (computeNMax(v.size()) != nVector->getNMax())
      nVector->setNMax(computeNMax(v.size()));
    NumVector<T>::operator=(v);
  }

  /// Copy operator.
  void operator= (const CoefficientVector<T>& cv) {
    if (cv.nVector->getNMax() != nVector->getNMax())
      nVector->setNMax(cv.nVector->getNMax());
    NumVector<T>::operator=(cv.getNumVector());
  } 

  const T& operator()(unsigned int i) const {
    return NumVector<T>::operator()(i);
  }
  T& operator()(unsigned int i) {
    return NumVector<T>::operator()(i);
  }

  /// const Matrix-like access operator.
  const T& operator()(int i, int j) const {
    return NumVector<T>::operator()(nVector->getIndex(i,j));
  }

  /// Matrix-like access operator.
  T& operator()(int i, int j) {
    return NumVector<T>::operator()(nVector->getIndex(i,j));
  }

  /// Set new coefficients in matrix form.
  /// If \f$n_{max}\f$ of new coefficient matrix is different, 
  /// CoefficientVector<T> will change accordingly.
  void setCoeffs(const NumMatrix<T>& coeffMatrix) {
    // set up IndexVector
    int nmax = coeffMatrix.getRows() - 1;
    if (nVector->getNMax() != nmax) {
      nVector->setNMax(nmax);
      // resize storage contained
      NumVector<T>::resize(nVector->getNCoeffs());
    }
    for (int n = 0; n < NumVector<T>::size(); n++)
      NumVector<T>::operator()(n) = coeffMatrix(nVector->getIndex1(n),nVector->getIndex2(n));
  }

  /// Map \p *this onto a matrix.
  void fillCoeffMatrix(NumMatrix<T>& coeffMatrix) const {
    // resize matrix if necessary
    int nmax = nVector->getNMax();
    if (coeffMatrix.getRows() != nmax + 1 || coeffMatrix.getColumns() != nmax +1)
      coeffMatrix.resize(nmax+1,nmax+1);
    coeffMatrix.clear();
    for (int n = 0; n < NumVector<T>::size(); n++)
      coeffMatrix(nVector->getIndex1(n),nVector->getIndex2(n)) = NumVector<T>::operator()(n);
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
    if (nVector->getNMax() != nmax) {
      nVector->setNMax(nmax);
      NumVector<T>::resize(nVector->getNCoeffs());
      NumVector<T>::clear();
    }
  }

  /// Get \f$n_{max}\f$.
  unsigned int getNMax() const {
    return nVector->getNMax();
  }
  /// Get \f$n_{coeffs}\f$.
  unsigned int getNCoeffs() const {
    return nVector->getNCoeffs();
  }

  /// Get the IndexVector vector <-> matrix mapping.
  const IndexVector& getIndexVector() const {
    return *nVector;
  }
  const NumVector<T>& getNumVector() const {
    return *this;
  }
 private:
  IndexVector* nVector;
  int computeNMax(int nCoeffs) {
    return (int)floor(-1.5 + 0.5*sqrt(9.0 + 8*(nCoeffs-1))); 
  }
  bool isComplex() {
    T ownval;
    complex<data_t> complexval;
    if (typeid(ownval) == typeid(complexval))
      return true;
    else
      return false;
  }
};

template <class T>
NumVector<T> operator*(const NumMatrix<T>& M, const CoefficientVector<T>& cv) {
  return M*cv.getNumVector();
}
#endif
