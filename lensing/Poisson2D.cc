// Function to invert tridiagonal matrix equation.
// Left, centre, and right diagonal elements of matrix
// stored in arrays a, b, c, respectively.
// Right-hand side stored in array w.
// Solution written to array u.
// Matrix is NxN. Arrays a, b, c, w, u assumed to be of extent N+2,
// with redundant 0 and N+1 elements.

#include <NumMatrix.h>
#include <NumVector.h>
#include <complex.h>

#ifndef FFT_CC
#include "FFT.cc"
typedef complex<double> Complex;
#endif
#define POISSON2D_CC

void Tridiagonal (NumVector<double>& a, NumVector<double>& b, NumVector<double>& c, NumVector<double>& w, NumVector<double>& u) {
  // Find N. Declare local arrays.
  int N = a.size() - 2;
  NumVector<double> x(N), y(N);
  // Scan up diagonal from i = N to 1
  x(N-1) = - a(N) / b(N);
  y(N-1) = w(N) / b(N);
  for (int i = N-2; i > 0; i--) {
    x(i) = - a(i+1) / (b(i+1) + c(i+1) * x(i+1));
    y(i) = (w(i+1) - c(i+1) * y(i+1)) / (b(i+1) + c(i+1) * x(i+1));
  }
  x(0) = 0;
  y(0) = (w(1) - c(1) * y(1)) / (b(1) + c(1) * x(1));
  // Scan down diagonal from i = 1 to N
  u(1) = y(0);
  for (int i = 1; i < N; i++)
    u(i+1) = x(i) * u(i) + y(i);
}

void Poisson2D (NumMatrix<double>& u, NumMatrix<double>& v, double alphaL, double betaL, NumVector<double>& gammaL, double alphaH, double betaH, NumVector<double>& gammaH, double dx, double L, int Neumann) {
  double kappa = M_PI * dx / L;
  // Find N and J. Declare local arrays.
  int N = u.getRows() - 2;
  int J = u.getColumns() - 1;
  NumMatrix<double> V(N+2, J+1), U(N+2, J+1);
  NumVector<double> GammaL(J+1), GammaH(J+1);
  // Fourier transform boundary conditions
  if (Neumann) {
    fft_cos_forward (gammaL, GammaL);
    fft_cos_forward (gammaH, GammaH);
  } else {
    fft_sin_forward (gammaL, GammaL);
    fft_sin_forward (gammaH, GammaH);
  }
  // Fourier transform source term
  NumVector<double> In(J+1);
  NumVector<double> Out(J+1);
  // FIXME: why starting with 1?
  for (int i = 1; i <= N; i++)
    {
      for (int j = 0; j <= J; j++) In(j) = v(i, j);
      if (Neumann)
	fft_cos_forward(In, Out);
      else
	fft_sin_forward(In, Out);
      for (int j = 0; j <= J; j++) V(i, j) = Out(j);
    }

  // Solve tridiagonal matrix equations
  if (Neumann)
    {
      for (int j = 0; j <= J; j++)
	{
	  NumVector<double> a(N+2), b(N+2), c(N+2), w(N+2), uu(N+2);
	  // Initialize tridiagonal matrix
	  for (int i = 2; i <= N; i++) a(i) = 1.;
	  for (int i = 1; i <= N; i++)
	    b(i) = -2. - (j*j) * kappa * kappa;
	  // modification: with alphaL = alphaH = 0
	  // U(i,j) = inf for j=0
	  // not shure if this is correct!
	  if (j != 0 || alphaL != 0)
	    b(1) -= betaL / (alphaL * dx - betaL);
	  // again modified
	  if (j != 0 || alphaH != 0)
	    b(N) += betaH / (alphaH * dx + betaH);
	  for (int i = 1; i <= N-1; i++) c(i) = 1.;
	  // Initialize right-hand side vector
	  for (int i = 1; i <= N; i++)
	    w(i) = V(i, j) * dx * dx;
	  w(1) -= GammaL(j) * dx / (alphaL * dx - betaL);
	  w(N) -= GammaH(j) * dx / (alphaH * dx + betaH);
	  // Invert tridiagonal matrix equation
	  Tridiagonal (a, b, c, w, uu);
	  for (int i = 1; i <= N; i++) {
	    U(i, j) = uu(i);
	    //if (j==0) std::cout << i << " " << uu(i) << "\t";
	  }
	}
    }
  else
    {
      for (int j = 1; j < J; j++)
	{
	  NumVector<double> a(N+2), b(N+2), c(N+2), w(N+2), uu(N+2);
	  // Initialize tridiagonal matrix
	  for (int i = 2; i <= N; i++) a(i) = 1.;
	  for (int i = 1; i <= N; i++)
	    b(i) = -2. - (j*j) * kappa * kappa;
	  b(1) -= betaL / (alphaL * dx - betaL);
	  b(N) += betaH / (alphaH * dx + betaH);
	  for (int i = 1; i <= N-1; i++) c(i) = 1.;
	  // Initialize right-hand side vector
	  for (int i = 1; i <= N; i++)
	    w(i) = V(i, j) * dx * dx;
	  w(1) -= GammaL(j) * dx / (alphaL * dx - betaL);
	  w(N) -= GammaH(j) * dx / (alphaH * dx + betaH);
	  // Invert tridiagonal matrix equation
	  Tridiagonal (a, b, c, w, uu);
	  for (int i = 1; i <= N; i++) U(i, j) = uu(i);
	}
      for (int i = 1; i <= N ; i++)
	{
	  U(i, 0) = 0.; U(i, J) = 0.;
	}
    }
  // Reconstruct solution via inverse Fourier transform
  for (int i = 1; i <= N; i++)
    {
      for (int j = 0; j <= J; j++) Out(j) = U(i, j);
      if (Neumann) 
	fft_cos_inverse (Out,In);
      else
	fft_sin_inverse (Out,In);
      for (int j = 0; j <= J; j++) u(i, j) = In(j);
    }
  // Calculate i=0 and i=N+1 values
  for (int j = 0; j <= J; j++)
    {
      u(0, j) = (gammaL(j) * dx - betaL * u(1, j)) /
	(alphaL * dx - betaL);
      u(N+1, j) = (gammaH(j) * dx + betaH * u(N, j)) /
	(alphaH * dx + betaH);
    }
 }
