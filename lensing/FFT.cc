#include <NumVector.h>
#include <NumMatrix.h>
#include <complex.h>
#include <fftw3.h>
#define FFT_CC

#ifndef POISSON2D_CC

typedef complex<double> Complex;
const Complex I = Complex(0,1);

void fft_forward(NumVector<double>& f, NumVector<Complex>& F) {
  int N = f.size();
  // data set must have even number of points
  if (N%2 ==0) {
    double *in;
    fftw_complex *out;
    fftw_plan p;
    in = reinterpret_cast<double(*)>(fftw_malloc(sizeof(fftw_complex)*N));
    out = reinterpret_cast<double(*)[2]>(fftw_malloc(sizeof(fftw_complex)*(N/2+1)));
    // copy NumVector into double*
    for (int i = 0; i < N; i++)
      in[i] = f(i);

    p = fftw_plan_dft_r2c_1d(N,in,out,FFTW_ESTIMATE);
    fftw_execute(p);
    // copy fftw_complex into NumVector
    if (F.size() != N) F = NumVector<Complex>(N);
    for (int i = 0; i < N/2+1; i++)
      F(i) = out[i][0] + I*out[i][1];
    // reconstruct the rest of the array from the hermiticity condition
    for (int i = 1; i < N/2; i++)
      F(N/2+i) = out[N/2-i][0] - I*out[N/2-i][1];
    fftw_destroy_plan(p);
    fftw_free(in);  
    fftw_free(out);
  }
}

void fft_forward(NumVector<Complex>& f, NumVector<Complex>& F) {
  int N = f.size();
  fftw_complex *data;
  fftw_plan p;
  data = reinterpret_cast<double(*)[2]>(fftw_malloc(sizeof(fftw_complex) * N));
  // copy NumVector into fftw_complex
  for (int i = 0; i < N; i++) {
    data[i][0] = real(f(i));
    data[i][1] = imag(f(i));
  }
  p = fftw_plan_dft_1d(N, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  // copy fftw_complex into NumVector
  if (F.size() != N) F = NumVector<Complex>(N);
  for (int i = 0; i < N; i++)
    F(i) = data[i][0] + I*data[i][1];

  fftw_destroy_plan(p);
  fftw_free(data);
}

void fft_forward(NumMatrix<double>& f, NumMatrix<Complex>& F) {
  int N = f.getRows();
  int J = f.getColumns();
  // data set must have even number of points in the last dimension
  if (J%2 ==0) {
    double *in;
    fftw_complex *out;
    fftw_plan p;
    in = reinterpret_cast<double(*)>(fftw_malloc(sizeof(fftw_complex)*N*J));
    out = reinterpret_cast<double(*)[2]>(fftw_malloc(sizeof(fftw_complex)*N*(J/2+1)));
    //copy NumVector into double*
    for (int i = 0; i < N; i++)
      for (int j = 0; j < J; j++)
	in[j+J*i] = f(i,j);

    p = fftw_plan_dft_r2c_2d(N,J,in,out,FFTW_ESTIMATE);
    fftw_execute(p);
    // copy fftw_complex into NumMatrix
    if (F.getRows() != N || F.getColumns() != J) F = NumMatrix<Complex>(N,J);
    int i_;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < J/2+1; j++)
	F(i,j) = out[j+(J/2+1)*i][0] + I*out[j+(J/2+1)*i][1];
      // reconstruct the rest of the array from the hermiticity condition
      for (int j = 1; j < J/2; j++) {
	if (i>0) i_ = N-i;
	else i_ = 0;
	F(i,J/2+j) = out[J/2-j+(J/2+1)*i_][0] - I*out[J/2-j+(J/2+1)*i_][1];
      }
    }
    fftw_destroy_plan(p);
    fftw_free(in);  
    fftw_free(out);
  }
}
    
void fft_forward(NumMatrix<Complex>& f, NumMatrix<Complex>& F) {
  int N = f.getRows();
  int J = f.getColumns();
  fftw_complex *data;
  fftw_plan p;
  data = reinterpret_cast<double(*)[2]>(fftw_malloc(sizeof(fftw_complex) * N*J));
  // copy NumMatrix into fftw_complex
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < J; j++) {
      data[j+J*i][0] = real(f(i,j));
      data[j+J*i][1] = imag(f(i,j));
    }
  }
  p = fftw_plan_dft_2d(N, J, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  // copy fftw_complex into NumVector
  if (F.getRows() != N || F.getColumns() != J) F = NumMatrix<Complex>(N,J);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < J; j++)
      F(i,j) = data[j+J*i][0] + I*data[j+J*i][1];

  fftw_destroy_plan(p);
  fftw_free(data);
}

void fft_sin_forward(NumVector<double>& f, NumVector<double>& F) {
  // because of the strange definition of the DST in fftw
  // we have to consrtuct the DST from the DFT
  int J = f.size() - 1;
  int N = 2*J;
  NumVector<Complex> Fcomplex(N), fcomplex(N);
  for (int i =1; i < J; i++) {
    fcomplex(i) = f(i);
    fcomplex(N-i) = -f(i);
  }
  fft_forward(fcomplex,Fcomplex);
  if (F.size() != f.size())
    F = NumVector<double>(f.size());
  for (int i =1; i < J; i++)
    F(i) = -2 * imag(Fcomplex(i));
  // don't know why but we need factor 4 here
  F /= 4;
}

void fft_cos_forward(NumVector<double>& f, NumVector<double>& F) {
//   int N = f.size();
//   double *data;
//   fftw_plan p;
//   data = reinterpret_cast<double(*)>(fftw_malloc(sizeof(fftw_complex) * N));
//   // copy NumVector into double*
//   for (int i = 0; i < N; i++)
//     data[i] = f(i);

//   p = fftw_plan_r2r_1d(N, data, data,FFTW_REDFT00, FFTW_ESTIMATE);
//   fftw_execute(p);
//   // copy double* into NumVector
//   if (F.size() != N) F = NumVector<double>(N);
//   for (int i = 0; i < N; i++)
//     F(i) = data[i];
//   fftw_destroy_plan(p);
//   fftw_free(data);  
  int J = f.size() - 1;
  int N = 2*J;
  NumVector<Complex> Fcomplex(N), fcomplex(N);
  fcomplex(0) = f(0);
  fcomplex(J) = f(J);
  for (int i =1; i < J; i++) {
    fcomplex(i) = fcomplex(N-i) = f(i);
  }
  fft_forward(fcomplex,Fcomplex);
  if (F.size() != f.size())
    F = NumVector<double>(f.size());
  F(0) = real(Fcomplex(0));
  F(J) = real(Fcomplex(J));
  for (int i =1; i < J; i++)
    F(i) = 2*real(Fcomplex(i));
 }

void fft_inverse(NumVector<Complex>& F, NumVector<double>& f) {
  int N = F.size();
  if (N%2 == 0) {
    double *out;
    fftw_complex *in;
    fftw_plan p;
    out = reinterpret_cast<double(*)>(fftw_malloc(sizeof(fftw_complex)*N));
    in = reinterpret_cast<double(*)[2]>(fftw_malloc(sizeof(fftw_complex)*N));
    // copy NumVector into fftw_complex*
    for (int i = 0; i < N/2+1; i++) {
      in[i][0] = real(F(i));
      in[i][1] = imag(F(i));
    }

    p = fftw_plan_dft_c2r_1d(N,in,out,FFTW_ESTIMATE);
    fftw_execute(p);
    // copy fftw_complex into NumVector
    if (f.size() != N) f = NumVector<double>(N);
    for (int i = 0; i < N; i++)
      f(i) = out[i]/N;

    fftw_destroy_plan(p);
    fftw_free(in);  
    fftw_free(out);
  }
}

void fft_inverse(NumVector<Complex>& F, NumVector<Complex>& f) {
  int N = F.size();
  fftw_complex *data;
  fftw_plan p;
  data = reinterpret_cast<double(*)[2]>(fftw_malloc(sizeof(fftw_complex) * N));
  // copy NumVector into fftw_complex
  for (int i = 0; i < N; i++) {
    data[i][0] = real(F(i));
    data[i][1] = imag(F(i));
  }
  p = fftw_plan_dft_1d(N, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  // copy fftw_complex into NumVector
  if (f.size() != N) f = NumVector<Complex>(N);
  for (int i = 0; i < N; i++)
    f(i) = 1./N *(data[i][0] + I*data[i][1]);

  fftw_destroy_plan(p);
  fftw_free(data);
}

void fft_inverse(NumMatrix<Complex>& F, NumMatrix<double>& f) {
  int N = F.getRows();
  int J = F.getColumns();
  if (J%2 == 0) {
    double *out;
    fftw_complex *in;
    fftw_plan p;
    out = reinterpret_cast<double(*)>(fftw_malloc(sizeof(fftw_complex)*N*J));
    in = reinterpret_cast<double(*)[2]>(fftw_malloc(sizeof(fftw_complex)*N*(J/2+1)));
    // copy NumVector into fftw_complex*
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < J/2+1; j++) {
	in[j+(J/2+1)*i][0] = real(F(i,j));
	in[j+(J/2+1)*i][1] = imag(F(i,j));
      }
    }

    p = fftw_plan_dft_c2r_2d(N,J,in,out,FFTW_ESTIMATE);
    fftw_execute(p);
    // copy fftw_complex into NumVector
    if (f.getRows() != N || f.getColumns() != J) f = NumMatrix<double>(N,J);
    
    for (int i = 0; i < N; i++)
    for (int j = 0; j < J; j++)
      f(i,j) = out[j+J*i]/(N*J);

    fftw_destroy_plan(p);
    fftw_free(in);  
    fftw_free(out);
  }
}

void fft_inverse(NumMatrix<Complex>& F, NumMatrix<Complex>& f) {
  int N = F.getRows();
  int J = F.getColumns();
  fftw_complex *in;
  fftw_complex *out;
  fftw_plan p;
  in = reinterpret_cast<double(*)[2]>(fftw_malloc(sizeof(fftw_complex) * N*J));
  out = reinterpret_cast<double(*)[2]>(fftw_malloc(sizeof(fftw_complex) * N*J));
  // copy NumMatrix into fftw_complex
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < J; j++) {
      in[j+J*i][0] = real(F(i,j));
      in[j+J*i][1] = imag(F(i,j));
    }
  }
  p = fftw_plan_dft_2d(N, J, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  // copy fftw_complex into NumVector
  if (f.getRows() != N || f.getColumns() != J) f = NumMatrix<Complex>(N,J);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < J; j++)
      f(i,j) = 1./(N*J) * (out[j+J*i][0] + I*out[j+J*i][1]);

  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
}

void fft_sin_inverse(NumVector<double>& F, NumVector<double>& f) {
  // same as the DST: we have to go over the DFT^(-1) to get
  // the DST^(-1)
  int J = F.size() - 1;
  int N = 2*J;
  NumVector<Complex> Fcomplex(N);
  NumVector<double> tmp(N);
  for (int i = 1; i < J; i++) {
    Fcomplex(i) = -2.*I*F(i);
    Fcomplex(N-i) = 2.*I*F(i);
  }
  fft_inverse(Fcomplex,tmp);
  if (f.size() != F.size()) f = NumVector<double>(F.size());
  f(0) = f(J) = 0;
  for (int i =1; i < J; i++) 
    f(i) = tmp(i);
}

void fft_cos_inverse(NumVector<double>& F, NumVector<double>& f) {
  int J = F.size() - 1;
  int N = 2*J;
  NumVector<Complex> Fcomplex(N);
  NumVector<double> tmp(N);
  Fcomplex(0) = F(0);
  Fcomplex(J) = F(J);
  for (int i = 1; i < J; i++) {
    Fcomplex(i) = Fcomplex(N-i) = F(i)/2;
  }
  fft_inverse(Fcomplex,tmp);
  if (f.size() != F.size()) f = NumVector<double>(F.size());
  for (int i =0; i <= J; i++) 
    f(i) = tmp(i);
}

double wavenumber(int j, double dx, int N) {
  if (j >= 0 && j < N/2) return 2*M_PI*j/(N*dx);
  if (j==N/2) return 2*M_PI/(2*dx);
  if (j>N/2 && j < N) return 2*M_PI*(j-N)/(N*dx);
//   if (j >= 0 && j < N/2) return 1.*j/(N*dx);
//   if (j==N/2) return 1./(2*dx);
//   if (j>N/2 && j < N) return 1.*(j-N)/(N*dx);
  else return 0;
}

#endif
