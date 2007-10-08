#ifndef PIXELCOVARIANCEMATRIX_H
#define PIXELCOVARIANCEMATRIX_H

#include <NumVector.h>
#include <NumMatrix.h>
#include <Typedef.h>
#include <History.h>
#include <frame/SegmentationMap.h>
#include <frame/CorrelationFunction.h>

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
  void setCovarianceMatrix(const CorrelationFunction& xi, const Grid& grid, History& history);
  /// Get the bandwidth of the covariance matrix.
  /// <tt>bandwidth</tt> is the number of bands in the matrix with values different from 0.
  unsigned int getBandwidth() const ;
  /// Set the bandwidth of the matrix.
  /// Existing entries are cleared after calling the method.
  void setBandwidth(unsigned int bandwidth);
  /// Set the entry of the given <tt>band</tt>.
  /// This set the values of the band with a given horizontal <tt>offset</tt> to
  /// <tt>val</tt>.
  void setBand(unsigned int band, int offset, data_t val);
  /// Get the offset of the given <tt>band</tt>.
  int getOffset(unsigned int band) const;
  /// Get the value of the given <tt>band</tt>.
  data_t getValue(unsigned int band) const;
  /// Index operator.
  /// Since the matrix is not stored as a whole but rather only the values and offsets 
  /// of all bands, this method returns <tt>val</tt> if <tt>i-j=offset</tt>.
  data_t operator()(unsigned int i, unsigned int j) const;
  /// Get a block matrix representation of dimension \f$N\times N\f$.
  NumMatrix<data_t> getMatrix(unsigned int N) const;
  /// Multiply with block matrix.
  /// The PixelCovarianceMatrix will have the appropriate number of columns
  /// to allow a defined multiplication. Thus, the resulting matrix will have the 
  /// dimensions of the matrix <tt>M</tt>.
  NumMatrix<data_t> operator*(const NumMatrix<data_t>& M) const;
  /// Multiply with NumVector.
  /// The PixelCovarianceMatrix will have the appropriate number of columns
  /// to allow a defined multiplication. Thus, the resulting vector will have the 
  /// size of the vector <tt>v</tt>.
  NumVector<data_t> operator*(const NumVector<data_t>& v) const;
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
  /// Save covariance matrix to text file.
  /// The format is as follows:
  /// \code
  /// # int band int offset data_t entry
  /// 0 -1 0.256
  /// 1  0 1.001
  /// 2  1 0.310
  /// \endcode
  void save(std::string filename) const;
  /// Load covariance matrix from a file created by save().
  void load(std::string filename);
    
 private:
  unsigned int bandwidth;
  NumVector<int> offset;
  NumVector<data_t> entry;
};

NumMatrix<data_t> operator*(const NumMatrix<data_t>& M, const PixelCovarianceMatrix& V);

#endif
