#ifndef PIXELCOVARIANCEMATRIX_H
#define PIXELCOVARIANCEMATRIX_H

#include <NumVector.h>
#include <NumMatrix.h>
#include <SegmentationMap.h>
#include <History.h>

/// Class for the pixel covariance matrix.
/// 
/// This class provides means to effectively store and work with pixel covariance matrices.
/// It is assumed here that the matrices are banded and the entries on each band are 
/// constant (which is true in the case of <tt>GAUSSIAN</tt> noise).
/// 
/// \todo
/// - implement non-constant entries
/// - exploit symmetry if present

class PixelCovarianceMatrix {
 public:
  /// Default constructor.
  PixelCovarianceMatrix();
  /// Argumented constructor.
  PixelCovarianceMatrix(unsigned int bandwidth);
  /// Argumented constuctor for estimating the covariance matrix from an Object.
  /// The constructor calls getCorrelationFunction() to measure the correlation function
  /// in dependence of the distance \f$\xi(r)\f$ in the pixel data.
  /// It's assumed that \f$\xi(r\geq2)\approx 0\f$, which is true for most images, 
  /// such that the maximum bandwidth is limited to 9 (the pixel itself and his 
  /// 8 neighbors).
  void setCovarianceMatrix(const NumVector<double>& data, const SegmentationMap& segMap, History& history);
  /// Get the bandwidth of the covariance matrix.
  /// <tt>bandwidth</tt> is the number of bands in the matrix with values different from 0.
  unsigned int getBandwidth() const ;
  /// Set the bandwidth of the matrix.
  /// Existing entries are cleared after calling the method.
  void setBandwidth(unsigned int bandwidth);
  /// Set the entry of the given <tt>band</tt>.
  /// This set the values of the band with a given horizontal <tt>offset</tt> to
  /// <tt>val</tt>.
  void setBand(unsigned int band, int offset, double val);
  /// Index operator.
  /// Since the matrix is not stored as a whole but rather only the values and offsets 
  /// of all bands, this method returns <tt>val</tt> if <tt>i-j=offset</tt>.
  double operator()(unsigned int i, unsigned int j) const;
  /// Get a block matrix representation of dimension \f$N\times N\f$.
  NumMatrix<double> getMatrix(unsigned int N) const;
  /// Multiply with block matrix.
  /// The PixelCovarianceMatrix will have the appropriate number of columns
  /// to allow a defined multiplication. Thus, the resulting matrix will have the 
  /// dimensions of the matrix <tt>M</tt>.
  NumMatrix<double> operator*(const NumMatrix<double>& M) const;
  /// Multiply with NumVector.
  /// The PixelCovarianceMatrix will have the appropriate number of columns
  /// to allow a defined multiplication. Thus, the resulting vector will have the 
  /// size of the vector <tt>v</tt>.
  NumVector<double> operator*(const NumVector<double>& v) const;
  /// Invert the covariance matrix.
  /// Since inversion of banded matrices is not easily done analytically, 
  /// a small portion of the covariance matrix is build as block matrix.
  /// The block matrix is then inverted numerically and the resulting offsets and 
  /// values from the inverse build up the new PixelCovarianceMatrix.\n
  /// Here entries with a smaller absolute value than <tt>0.001*diag</tt>, 
  /// where <tt>diag</tt> is the value on the diagonal, are neglected to keep the
  /// number of bands of the inverse reasonably small.
  PixelCovarianceMatrix invert() const;
  /// Transpose covariance matrix.
  PixelCovarianceMatrix transpose() const;
  /// Measures the correlation function <tt>corr</tt> in dependence of the 
  /// <tt>distance</tt> in the pixel data of the Object.
  /// It uses the SegmentationMap of the Object to decide which pixels belong to 
  /// the noise. Thus, this method provides the pixel correlation function of 
  /// all noise pixels in the image.
  void getCorrelationFunction(const NumVector<double>& data, const SegmentationMap& segMap, NumVector<double>& corr, NumVector<double>& distance) const;
 private:
  unsigned int bandwidth;
  NumVector<int> offset;
  NumVector<double> entry;
};

NumMatrix<double> operator*(const NumMatrix<double>& M, const PixelCovarianceMatrix& V);

#endif
