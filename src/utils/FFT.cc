#include "../../include/utils/FFT.h"
#include <stdexcept>

using namespace shapelens;
typedef complex<data_t> Complex;

FourierTransform1D::FourierTransform1D() : NumVector<Complex> () {
  N = 0;
} 
FourierTransform1D::FourierTransform1D(unsigned int N) : 
  NumVector<Complex> (N/2+1), N(N) {}

complex<data_t>& FourierTransform1D::operator()(unsigned int i) {
  if (i <= N/2)
    return NumVector<Complex>::operator()(i);
  else
    return NumVector<Complex>::operator()(N-(int)i);
    
}
const complex<data_t>& FourierTransform1D::operator()(unsigned int i) const {
  if (i <= N/2)
    return NumVector<Complex>::operator()(i);
  else
    return NumVector<Complex>::operator()(N-(int)i);
}
data_t FourierTransform1D::getWavenumber(int i) const {
  if (i >= 0 && i < N/2) return 2*M_PI*i/N;
  if (i==N/2) return M_PI;
  if (i>N/2 && i < N) return 2*M_PI*(i-N)/N;
}
void FourierTransform1D::resize(unsigned int N1) {
  N = N1;
  NumVector<Complex>::resize(N/2+1);
}

int FourierTransform1D::getRealSize() const {
  return N;
}

FourierTransform2D::FourierTransform2D() : NumVector<Complex> () {
  N = J = 0;
} 
FourierTransform2D::FourierTransform2D(unsigned int N, unsigned int J) : 
  NumVector<Complex> (N*(J/2+1)), N(N), J(J) {}

complex<data_t>& FourierTransform2D::operator()(unsigned int i, unsigned int j) {
  if (j <= J/2)
    return NumVector<Complex>::operator()(i*(J/2+1) + j);
  else
    return NumVector<Complex>::operator()(i*(J/2+1) + (J-(int)j));
}

const complex<data_t>& FourierTransform2D::operator()(unsigned int i, unsigned int j)  const {
  if (j <= J/2)
    return NumVector<Complex>::operator()(i*(J/2+1) + j);
  else
    return NumVector<Complex>::operator()(i*(J/2+1) + (J-(int)j));
}

int FourierTransform2D::getRealSize(bool dimension) const {
  if (!dimension)
    return N;
  else
    return J;
}

complex<data_t> FourierTransform2D::getWavenumber(int i, int j) const {
  return complex<data_t>(wavenumber(i,0),wavenumber(j,1));
}

data_t FourierTransform2D::wavenumber(int k, bool dimension) const {
  int K;
  if (dimension == false)
    K = N;
  else
    K = J;
  if (k >= 0 && k < K/2) return 2*M_PI*k/K;
  if (k==K/2) return M_PI;
  if (k>K/2 && k < K) return 2*M_PI*(k-K)/K;
}

void FourierTransform2D::resize(unsigned int N1, unsigned int J1) {
  N = N1;
  J = J1;
  NumVector<Complex>::resize(N*(J/2+1));
}


void FFT::transform(const NumVector<data_t>& f, FourierTransform1D& F) {
  int N = f.size();
  // data set must have even number of points
  if (N%2 ==0) {
    if (F.size() != N)
      F.resize(N);
    // we work directly on the contents of f and F
    fftw_plan p = fftw_plan_dft_r2c_1d(N,const_cast<data_t*>(f.c_array()),reinterpret_cast<data_t(*)[2]>(F.c_array()),FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  } else 
    throw std::length_error("FFT: NumVector has odd number of entries");
}

void FFT::transform(const FourierTransform1D& F, NumVector<data_t>& f) {
  int N = F.getRealSize();
  // data set must have even number of points
  if (N%2 ==0) {
    if (f.size() != N)
      f.resize(N);
    // we work directly on the contents of f and F
    // CAUTION: Due to Hermiticity of real -> Complex FFT, 
    // not all entries of F are filled
    fftw_plan p = fftw_plan_dft_c2r_1d(N,reinterpret_cast<data_t(*)[2]>(const_cast<Complex*>(F.c_array())),reinterpret_cast<data_t*>(f.c_array()),FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    f /= N;
  } else 
    throw std::length_error("FFT: FourierTransform1D has odd number of entries");
}

void FFT::transform(const NumMatrix<data_t>& f, FourierTransform2D& F) {
  int N = f.getRows();
  int J = f.getColumns();

  // data set must have even number of points in the last dimension
  if (J%2 ==0) {
    if (F.getRealSize(0) != N || F.getRealSize(1) != J)
      F.resize(N,J);
    // we work directly on the contents of f and F
    fftw_plan p = fftw_plan_dft_r2c_2d(N,J,const_cast<data_t*>(f.c_array()),reinterpret_cast<data_t(*)[2]>(F.c_array()),FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  } else 
    throw std::length_error("FFT: NumMatrix has odd number of columns");
}

void FFT::transform(const FourierTransform2D& F, NumMatrix<data_t>& f) {
  int N = F.getRealSize(0);
  int J = F.getRealSize(1);
  if (J%2 == 0) {
    if (f.getRows() != N || f.getColumns() != J)
      f.resize(N,J);
    // we work directly on the contents of f and F
    fftw_plan p = fftw_plan_dft_c2r_2d(N,J,reinterpret_cast<data_t(*)[2]>(const_cast<Complex*>(F.c_array())),f.c_array(),FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    // normalization
    f /= N*J;
  } else 
    throw std::length_error("FFT: FourierTransform2D has odd number of columns");
}

void FFT::transform(const Image<data_t>& f, FourierTransform2D& F) {
  int N = f.getSize(0);
  int J = f.getSize(1);

  // data set must have even number of points in the last dimension
  if (J%2 ==0) {
    if (F.getRealSize(0) != N || F.getRealSize(1) != J)
      F.resize(N,J);
    // we work directly on the contents of f and F
    fftw_plan p = fftw_plan_dft_r2c_2d(N,J,const_cast<data_t*>(f.c_array()),reinterpret_cast<data_t(*)[2]>(F.c_array()),FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  } else 
    throw std::length_error("FFT: Image has odd number of columns");
  
}

void FFT::transform(const FourierTransform2D& F, Image<data_t>& f) {
  int N = F.getRealSize(0);
  int J = F.getRealSize(1);
  if (J%2 == 0) {
    if (f.getSize(0) != N || f.getSize(1) != J) {
      f.resize(N,J);
      f.grid = Grid(0,0,N,J);
    }
    // we work directly on the contents of f and F
    fftw_plan p = fftw_plan_dft_c2r_2d(N,J,reinterpret_cast<data_t(*)[2]>(const_cast<Complex*>(F.c_array())),f.c_array(),FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    // normalization
    f /= N*J;
  } else 
    throw std::length_error("FFT: FourierTransform2D has odd number of columns");
}

void FFT::convolve(Image<data_t>& data, const Image<data_t>& kernel) {
  int M = data.getSize(0);
  int N = data.getSize(1);

  // FFT the data
  FourierTransform2D data_transformed(M,N);
  FFT::transform(data, data_transformed);
  FFT::convolve(data,data_transformed,kernel);
}

void FFT::convolve(Image<data_t>& data, FourierTransform2D& data_transformed, const Image<data_t>& kernel) {
  int M = data.getSize(0);
  int N = data.getSize(1);

  FourierTransform2D kernel_transformed(M,N);
  // ...and the convolution kernel
  // check size of kernel
  // if is doesn't fit to data: resize it
  if (N != kernel.getSize(0) || M != kernel.getSize(1)) {
    int N1 = kernel.getSize(0);
    int M1 = kernel.getSize(1); 
    if (N1%2==1)
      N1++;
    if (M1%2==1)
      M1++;
    int xmin = (N1-N)/2, xmax = N1 + (N-N1)/2, ymin = (M1-M)/2, ymax = M1+ (M-M1)/2;
    Image<data_t> kernel_resized(xmax-xmin,ymax-ymin);
    kernel.slice(kernel_resized,Point<int>(xmin,ymin),Point<int>(xmax,ymax));
    FFT::transform(kernel_resized, kernel_transformed);
  } 
  else
    FFT::transform(kernel, kernel_transformed);
    
  // make use of the convolution theorem in Fourier space:
  FFT::conv_multiply(data_transformed,kernel_transformed,data_transformed);
	
  // Transform back to real space
  FFT::transform(data_transformed,data);
}

void FFT::conv_multiply(const FourierTransform2D& f1, const FourierTransform2D& f2, FourierTransform2D& target) {
  int M = f1.getRealSize(0);
  int N = f1.getRealSize(1);
  for (int i=0; i<M; i++) {
    for (int j=0; j<N/2+1; j++) {
      target(i,j) = f1(i,j)*f2(i,j);
      // reorder the frequencies to have centered FFT afterwards
      if ((i+j)%2 == 1)
	target(i,j) *= -1;
    }
  }
}
