#include <utils/FFT.h>

#ifdef HAS_FFTW3
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
    return NumVector<Complex>::operator()(N/2-i);
    
}
const complex<data_t>& FourierTransform1D::operator()(unsigned int i) const {
  if (i <= N/2)
    return NumVector<Complex>::operator()(i);
  else
    return NumVector<Complex>::operator()(N/2-i);
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

FourierTransform2D::FourierTransform2D() : NumMatrix<Complex> () {
  N = J = 0;
} 
FourierTransform2D::FourierTransform2D(unsigned int N, unsigned int J) : 
  NumMatrix<Complex> (N,J/2+1),N(N),J(J) {}

complex<data_t>& FourierTransform2D::operator()(unsigned int i, unsigned int j) {
  if (j <= N/2)
    return NumMatrix<Complex>::operator()(i,j);
  else
    return NumMatrix<Complex>::operator()(i,N/2-j);
    
}

const complex<data_t>& FourierTransform2D::operator()(unsigned int i, unsigned int j)  const {
  if (j <= N/2)
    return NumMatrix<Complex>::operator()(i,j);
  else
    return NumMatrix<Complex>::operator()(i,N/2-j);
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
  NumMatrix<Complex>::resize(N,J/2+1);
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
  }
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
  }
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
  }
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
  }
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
  }
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
  }
}


#endif
