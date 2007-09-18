#ifndef POLARTRANSFORMATION_H
#define POLARTRANSFORMATION_H

#include <gsl/gsl_sf.h>
#include <NumVector.h>
#include <NumMatrix.h>
#include <Typedef.h>
#include <frame/Grid.h>
#include <shapelets/IndexVector.h>

/// Class for Polar Shapelets.
/// Performs a coefficient transformation cartesian -> polar and vice versa.
/// Therefore rearranges the coeffs matrices into vectors and transform via 
/// \f$w = A\cdot v\f$.
/// The transformation matrix and its inverse are stored since they only depend on
/// the maximal shapelet order \f$n_{max}\f$ and can therefore be reused.
/// 
/// See Paper III, section 2.3, for details.

class PolarTransformation {
 public:
  /// Default constructor.
  PolarTransformation();
  /// Argumented contructor.
  /// Transforms cartesian coefficients up to the order \f$n_1 + n_2 <= n_{max}\f$
  /// and polar coefficients up to order \f$n <= n_{max}\f$.
  PolarTransformation(unsigned int nmax);
  
  /// Return highest transformed order \f$n_{max}\f$.
  unsigned int getOrder();
  /// set the highest transformed order \f$n_{max}\f$.
  void setOrder (unsigned int nmax);
  /// Return polar coeffs for external transformations.
  void getPolarCoeffs(const NumMatrix<data_t>& cartesianCoeffs, NumMatrix<complex<data_t> >& polarCoeffs);
  /// Return cartesian coeffs for external transformations.
  void getCartesianCoeffs(const NumMatrix<complex<data_t> >& polarCoeffs, NumMatrix<data_t>& cartesianCoeffs);
 private:
  unsigned int nmax;
  NumMatrix<complex<data_t> > c2p,p2c;
  IndexVector nVector;
  void buildTransformationMatrix();
};

#endif
