#ifndef POLARTRANSFORMATION_H
#define POLARTRANSFORMATION_H

#include <gsl/gsl_sf.h>
#include <NumVector.h>
#include <NumMatrix.h>
#include <Typedef.h>
#include <frame/Grid.h>
#include <shapelets/CoefficientVector.h>

namespace shapelens {

/// Class for Polar Shapelets.
/// Performs a coefficient transformation cartesian -> polar and vice versa.
/// The transformation matrix and its inverse are stored since they only depend on
/// the maximal shapelet order \f$n_{max}\f$ and can therefore be reused.
/// 
/// See Paper III, section 2.3, for details.

class PolarTransformation {
 public:
  /// Default constructor.
  PolarTransformation();
  /// Compute \p polarCoeffs from \p cartesianCoeffs.
  /// The transformaition automatically adapts to the order of
  /// \p cartesianCoeffs.\n
  /// If \p covariance and \p polarCovariance are given, the transformation is also 
  /// applied to the covariance matrices.
  void getPolarCoeffs(const CoefficientVector<data_t>& cartesianCoeffs, CoefficientVector<complex<data_t> >& polarCoeffs, NumMatrix<data_t>* covariance = NULL, NumMatrix<complex<data_t> >* polarCovariance = NULL);
  /// Compute \p cartesianCoeffs from \p polarCoeffs.
  /// The transformaition automatically adapts to the order of
  /// \p polarCoeffs.
  /// If \p covariance and \p polarCovariance are given, the transformation is also 
  /// applied to the covariance matrices.
  void getCartesianCoeffs(const CoefficientVector<complex<data_t> >& polarCoeffs, CoefficientVector<data_t>& cartesianCoeffs, NumMatrix<complex<data_t> >* polarCovariance = NULL, NumMatrix<data_t>* covariance = NULL);
 private:
  unsigned int nmax;
  NumMatrix<complex<data_t> > c2p,p2c;
  void buildTransformationMatrix(const IndexVector& , const IndexVector& );
};
} // end namespace
#endif
