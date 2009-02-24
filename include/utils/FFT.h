#ifndef FFT_H
#define FFT_H

#ifdef HAS_FFTW3

#include <Typedef.h>
#include <fftw3.h>
#include <NumVector.h>
#include <NumMatrix.h>
#include <frame/Image.h>

/// Class for storing result of one-dimensional FFT.
/// This class exploits that the Fourier transform of a real vector
/// obeys the Hermiticity condition
/// \f[F(i) = F^*(N-i)\f]
/// for a vector index \f$i\f$ and the length of the real vector \f$N\f$.
class FourierTransform1D : public NumVector<complex<data_t> > {
 public:
  /// Constructor.
  FourierTransform1D();
  /// Constructor for real vector of size \f$N\f$.
  FourierTransform1D(unsigned int N);
  /// Access operator.
  complex<data_t>& operator()(unsigned int i);
  /// Access operator.
  const complex<data_t>& operator()(unsigned int i) const;
  /// Get wavenumber of index \p i of the transform.
  data_t getWavenumber(int i) const;
  /// Resize transform for real vector of size \f$N\f$.
  void resize(unsigned int N);
  /// Get size of real vector.
  int getRealSize() const;
 private:
  int N;
};

/// Class for storing result of one-dimensional FFT.
/// This class exploits that the Fourier transform of a real matrix
/// obeys the Hermiticity condition
/// \f[F(i,j) = F^*(i,J-j)\f]
/// for vector indices \f$i, j\f$ and the column number  of the real matrix \f$J\f$.
class FourierTransform2D : public NumMatrix<complex<data_t> > {
 public:
  /// Constructor.
  FourierTransform2D();
  /// Constructor for real matrix of size \f$N\times J\f$.
  FourierTransform2D(unsigned int N, unsigned int J);
  /// Access operator.
  complex<data_t>& operator()(unsigned int i, unsigned int j);
  /// Access operator.
  const complex<data_t>& operator()(unsigned int i, unsigned int j) const;
  /// Get wavenumber of index \p i of the transform.
  complex<data_t> getWavenumber(int i, int j) const;
  /// Resize transform for real matrix of size \f$N\times J\f$.
  void resize(unsigned int N, unsigned int J);
  /// Get size of real matrix in \p dimension.
  int getRealSize(bool dimension) const;
 private:
  int N,J;
  data_t wavenumber(int k, bool dimension) const;
};

/// class for one- and two-dimensional Fourier Transforms.
class FFT {
 public:
  /// Transform \f$f(x)\rightarrow F(k)\f$.
  static void transform(const NumVector<data_t>& f, FourierTransform1D& F);
  /// Transform \f$F(k)\rightarrow f(x)\f$.
  static void transform(const  FourierTransform1D& F, NumVector<data_t>& f);
  /// Transform two-dimensional \f$f(\vec{x})\rightarrow F(\vec{k})\f$.
  static void transform(const NumMatrix<data_t>& f,  FourierTransform2D& F);
  /// Transform two-dimensional \f$F(\vec{k})\rightarrow f(\vec{x})\f$.
  static void transform(const  FourierTransform2D& F, NumMatrix<data_t>& f);
  /// Transform Image \f$f(\vec{x})\rightarrow F(\vec{k})\f$.
  static void transform(const Image<data_t>& f,  FourierTransform2D& F);
  /// Transform Image \f$F(\vec{k})\rightarrow f(\vec{x})\f$.
  static void transform(const  FourierTransform2D& F, Image<data_t>& f);
  /// Convolve \p data with kernel
  /// It is assumed that both \p data and \p kernel have the same sizes.
  static void convolve(Image<data_t>& data, const Image<data_t>& kernel);
 private:
  friend class Object;
  static void reorder(Image<data_t>& im);
  static void conv_multiply(const FourierTransform2D& f1, const FourierTransform2D& f2, FourierTransform2D& target);
};

#endif // HAS_FFTW3
#endif // FFT_H
