#include "Poisson2D.cc"
#include "FFT.cc"
#include <gsl/gsl_math.h>

// compute potential for a given convergence map
// by solving poissons eq. d^2 psi = kappa with kappa = 0 on boundary
void computePotentialFromConvergence(NumMatrix<double>& kappa, NumMatrix<double>& psi, double dx, double L) {
  int N = kappa.getRows();
  int J = kappa.getColumns();
  if (psi.getRows() != N || psi.getColumns() != J) psi = NumMatrix<double>(N,J);
    // boundary conditions: potential at boundary is zero
  NumVector<double> gammaL(J),gammaH(J); 
  // solve poisson eq for Dirichlet boundary condition
  Poisson2D(psi,kappa,1,0,gammaL,1,0,gammaH,dx,L,0);
}

// compute shear from convergence with discrete fourier transforms and derivatives
// in fourier space
void computeShearFourier(NumMatrix<Complex>& Kappa, NumMatrix<Complex>& shear, double dx) {
  int N = Kappa.getRows();
  int J = Kappa.getColumns();
  if (shear.getRows() != N || shear.getColumns() != J)
    shear = NumMatrix<Complex>(N,J);
  NumMatrix<Complex> Shear1(N,J), Shear2(N,J);
  NumMatrix<double> tmp1(N,J), tmp2(N,J);
  double k1,k2;
  // 1. component
  for (int i=0; i<N; i++) {
    for (int j=0; j<J; j++) {
      k1 = wavenumber(i,dx,N);
      k2 = wavenumber(j,dx,J);
      if (i!=0 || j!=0) {
	Shear1(i,j) = (k1*k1-k2*k2)/(k1*k1+k2*k2) * Kappa(i,j);
	Shear2(i,j) = k1*k2/(2*(k1*k1+k2*k2)) * Kappa(i,j);
      }
    }
  }
  fft_inverse(Shear1,tmp1);
  fft_inverse(Shear2,tmp2);

  for(int i=0; i<N; i++)
    for (int j=0; j<J; j++)
      shear(i,j) = tmp1(i,j) + I*tmp2(i,j);
}

// compute shear from convergence with discrete fourier transforms and derivatives
// in fourier space
void computeFlexion1Fourier(NumMatrix<Complex>& Kappa, NumMatrix<Complex>& F, double dx) {
  int N = Kappa.getRows();
  int J = Kappa.getColumns();
  if (F.getRows() != N || F.getColumns() != J)
    F = NumMatrix<Complex>(N,J);
  NumMatrix<Complex> F1(N,J), F2(N,J);
  NumMatrix<double> tmp1(N,J), tmp2(N,J);
  double k1,k2;
  // 1. component
  for (int i=0; i<N; i++) {
    for (int j=0; j<J; j++) {
      k1 = wavenumber(i,dx,N);
      k2 = wavenumber(j,dx,J);
      if (i!=0 || j!=0) {
	F1(i,j) = I*k1 * Kappa(i,j);
	F2(i,j) = I*k2 * Kappa(i,j);
      }
    }
  }
  fft_inverse(F1,tmp1);
  fft_inverse(F2,tmp2);

  for(int i=0; i<N; i++)
    for (int j=0; j<J; j++)
      F(i,j) = tmp1(i,j) + I*tmp2(i,j);
}

void computeFlexion2Fourier(NumMatrix<Complex>& Kappa, NumMatrix<complex<double> >& G, double dx) {
  int N = G.getRows();
  int J = G.getColumns();
  if (G.getRows() != N || G.getColumns() != J)
    G = NumMatrix<Complex>(N,J);
  NumMatrix<complex<double> > G1(N,J), G2(N,J);
  NumMatrix<complex<double> > tmp1(N,J), tmp2(N,J);
  double k1,k2;

  for (int i=0;i < N; i++) {
    for(int j=0; j < J; j++) {
      k1 = wavenumber(i,dx,N);
      k2 = wavenumber(j,dx,J);
      if (i!=0 || j !=0) {
	G1(i,j) = I*(k1*k1*k1 - 3*k1*k2*k2)/(k1*k1+k2*k2) * Kappa(i,j);
	G2(i,j) = I*(3*k1*k1*k2 - k2*k2*k2)/(k1*k1+k2*k2) * Kappa(i,j);
      }
    }
  }
  fft_inverse(G1,tmp1);
  fft_inverse(G2,tmp2);
  for (int i=0;i < N; i++)
    for(int j=0; j < J; j++)
      G(i,j) = tmp1(i,j) + I*tmp2(i,j);
}

// compute convergence from pontential with finite differences
void computeConvergence(NumMatrix<double>& psi, NumMatrix<Complex>& kappa, double dx) {
  int J = psi.getColumns();
  int N = psi.getRows();
  if (kappa.getRows() != N || kappa.getColumns() != J)
    kappa = NumMatrix<Complex> (N,J);
  double x,y;
  for (int i = 1; i< N-1; i++) {
    for (int j = 1; j< N-1; j++) {
      x = -10 +i*dx;
      y = -10 +j*dx;
      kappa(i,j) =  (psi(i+1,j) - 2*psi(i,j) + psi(i-1,j))/2;
      kappa(i,j) += (psi(i,j+1) - 2*psi(i,j) + psi(i,j-1))/2;
      kappa(i,j) /= dx*dx/2;
    }
  }
}

// compute shear from potential with finite differences
void computeShear(NumMatrix<double>& psi, NumMatrix<Complex>& shear, double dx) {
  int J = psi.getColumns();
  int N = psi.getRows();
  if (shear.getRows() != N || shear.getColumns() != J)
    shear = NumMatrix<Complex> (N,J);
  double x,y;
  for (int i = 1; i< N-1; i++) {
    for (int j = 1; j< J-1; j++) {
      x = -10 +i*dx;
      y = -10 +j*dx;
      // 1. compontent
      shear(i,j) =  (psi(i+1,j) - 2*psi(i,j) + psi(i-1,j))/2;
      shear(i,j) -= (psi(i,j+1) - 2*psi(i,j) + psi(i,j-1))/2;
      // 2. component
      shear(i,j) += I*(psi(i+1,j+1) - psi(i+1,j-1) - psi(i-1,j+1) + psi(i-1,j-1))/4.;
      shear(i,j) /= dx*dx/2;
    }
  }
}

// compute 1. or 2. flexion from potential with finite differences
void computeFlexion(NumMatrix<double>& psi, bool number, NumMatrix<Complex>& flexion, double dx) {
  int J = psi.getColumns();
  int N = psi.getRows();
  if (flexion.getRows() != N || flexion.getColumns() != J)
    flexion = NumMatrix<Complex> (N,J);
  double psi111,psi222,psi112, psi122;
  for (int i = 3; i< N-3; i++) {
    for (int j = 3; j< J-3; j++) {
     // centered finite differences
      psi111 = 0.125*(psi(i+3,j) - 3*psi(i+1,j) + 3*psi(i-1,j) - psi(i-3,j));
      psi112 =  0.125*(psi(i+2,j+1)-2*psi(i,j+1) + psi(i-2,j+1));
      psi112 -= 0.125*(psi(i+2,j-1)-2*psi(i,j-1) + psi(i-2,j-1));
      psi122 = 0.125*(psi(i+1,j+2) - 2*psi(i+1,j) + psi(i+1,j-2));
      psi122 -= 0.125*(psi(i-1,j+2) - 2*psi(i-1,j) + psi(i-1,j-2));
      psi222 = 0.125*(psi(i,j+3) - 3*psi(i,j+1) + 3*psi(i,j-1) - psi(i,j-3));

      if (number==0) {
	// F1 = 1/2 (psi,111 + psi,122)
	// F2 = 1/2 (psi,112 + psi,222)
	flexion(i,j) = 0.5*(psi111 + psi122);
	flexion(i,j) += I*0.5*(psi112 + psi222);
      } else {
	// G1 = 1/2 (psi,111 - 3*psi,122)
	// G2 = 1/2 (3*psi,112 - psi,222)
	flexion(i,j) = 0.5*(psi111 - 3*psi122);
	flexion(i,j) += I*0.5*(3*psi112 - psi222);
      }
      flexion(i,j) /= dx*dx*dx/2;
    }
  }
}

// using Bartelmann/Scheider (2001), eg. 5.17
// d^2 kappa = grad F with n * grad kappa = n * F at the boundary
// we solve poisson eq. for grad F computed from derivatives of the shear
// by finite differences
// The von Neumann condition can be changed to a trivial Dirichlet condition
// as long as the convergence is 0 at the boundary
// If it is not, one has to use the von Neumann condition
// which does not work with this implementation of the poisson solver
// since it has du/dy = 0 at the left and the right boundary and arbitrary
// choices are not possible
// Use Successive Overrelaxation (Numerical recipes, p.857)
void computeConvergenceFromShear(NumMatrix<Complex>& shear, NumMatrix<Complex>& kappa, double L, double dx) {
  int N = shear.getRows();
  int J = shear.getColumns();
  if (kappa.getRows() != N || kappa.getColumns() != J)
    kappa = NumMatrix<Complex> (N,J);
  NumMatrix<double> rho(N,J), tmp(N,J);
  NumVector<double> gammaL(J), gammaH(J);
  // leaving out the first and the last index due to finite differences
  for (int i=1; i<N-1; i++) {
    for (int j=1; j<J-1; j++) {
      rho(i,j) =  real(shear(i+1,j))-2*real(shear(i,j))+real(shear(i-1,j));
      rho(i,j) -= real(shear(i,j+1))-2*real(shear(i,j))+real(shear(i,j-1));
      rho(i,j) += 0.5*(imag(shear(i+1,j+1)) - imag(shear(i+1,j-1)) - imag(shear(i-1,j+1)) + imag(shear(i-1,j-1)));
      rho(i,j) /= dx*dx;
    }
  }
  // von Neumann condition leeds to problems in case of noisy images
  // Poisson2D(tmp,Kappa,rho,0,1,gammaL,0,1,gammaH,dx,L,1);
  // using Dirichlet condition (see note above!)
  Poisson2D(tmp,rho,1,0,gammaL,1,0,gammaH,dx,L,0);
  for (int i=0; i<N; i++)
    for (int j=0; j<J; j++)
      kappa(i,j) = tmp(i,j);
}

// same as above but now compute grad F directly from F
// with finite differences
// this only works for first flexion
void computeConvergenceFromFlexion1(NumMatrix<Complex>& flexion,NumMatrix<Complex>& kappa, double L, double dx) {
  int N = flexion.getRows();
  int J = flexion.getColumns();
  if (kappa.getRows() != N || kappa.getColumns() != J)
    kappa = NumMatrix<Complex> (N,J);
  NumMatrix<double> rho(N,J), tmp(N,J);
  NumVector<double> gammaL(J), gammaH(J);
  // leaving out the the last index due to finite differences
  for (int i=1; i<N-1; i++) {
    for (int j=1; j<J-1; j++) {
      rho(i,j) = real(flexion(i+1,j)) - real(flexion(i-1,j));
      rho(i,j) += imag(flexion(i,j+1)) - imag(flexion(i,j-1));
      rho(i,j) /= dx*2;
    }
  }
  // same argument as above
  // Poisson2D(tmp,Kappa,rho,0,1,gammaL,0,1,gammaH,dx,L,1);
  Poisson2D(tmp,rho,1,0,gammaL,1,0,gammaH,dx,L,0);
  for (int i=0; i<N; i++)
    for (int j=0; j<J; j++)
      kappa(i,j) = tmp(i,j);
}


void computeConvergenceFromFlexionFourier(NumMatrix<Complex>& flexion, bool number, NumMatrix<Complex>& kappa, double dx) {
  int N = flexion.getRows();
  int J = flexion.getColumns();
  if (kappa.getRows() != N || kappa.getColumns() != J)
    kappa = NumMatrix<Complex> (N,J);
  NumMatrix<double> f1(N,J), f2(N,J), tmp(N,J);
  NumMatrix<Complex> Kappa(N,J), F1(N,J), F2(N,J), F(N,J);
  // should we exclude the borders for periodicity: borders = 0?
  for (int i=1; i<N-1; i++) {
    for(int j=1; j<J-1; j++) {
      f1(i,j) = real(flexion(i,j));
      f2(i,j) = imag(flexion(i,j));
    }
  }
  fft_forward(f1,F1);
  fft_forward(f2,F2);
  double k1,k2;
  for (int i=0; i<N; i++) {
    for(int j=0; j< J; j++) {
      k1 = wavenumber(i,dx,N);
      k2 = wavenumber(j,dx,J);
      if (number==0) {
	// from Bacon et al (2005)
	if (i!=0 || j!=0) Kappa(i,j) = -I/(k1*k1 + k2*k2) * (k1*F1(i,j) +k2*F2(i,j));
	// direct derivation using complex derivatives
	// has cross shaped artefacts because of exclusion of monopole?
	//if (i!=0) Kappa(i,j) += -I*F1(i,j)/k1;
	//if (j!=0) Kappa(i,j) += -I*F2(i,j)/k2;
      } else {
	// computed from complex derivatives
	//if (i!=0 || j!=0) {
	//  Kappa(i,j) =I*k1*F1(i,j) + k2*F1(i,j) - k1*F2(i,j) + I*k2*F2(i,j);
	//  Kappa(i,j)/= -k1*k1 -2.*I*k1*k2 +k2*k2;
	//}
	// naive solution in Fourier space
	//if (i!=0 && j!=0) {
	  // Kappa(i,j) = -I*(k1*k1+k2*k2)/(k1*k1*k1 - 3*k1*k2*k2)*F1(i,j);
	//  Kappa(i,j)+= -I*(k1*k1+k2*k2)/(3*k1*k1*k2 - k2*k2*k2)*F2(i,j);
	//}
	//with a minimizing the contribution at the divergences
	Kappa(i,j) = 0;
	if (i!=0) {  
	  double a = 1./(gsl_pow_2(3*k1*k1*k2-k2*k2*k2)/gsl_pow_2(k1*k1*k1-3*k1*k2*k2) + 1);
	  Kappa(i,j) += -a*I*(k1*k1+k2*k2)/(k1*k1*k1 - 3*k1*k2*k2)*F1(i,j);
	}
	if (j!=0) {
	  // b= 1-a
	  double b = 1./(gsl_pow_2(k1*k1*k1-3*k1*k2*k2)/gsl_pow_2(3*k1*k1*k2-k2*k2*k2) + 1);
	  Kappa(i,j)+= -b*I*(k1*k1+k2*k2)/(3*k1*k1*k2 - k2*k2*k2)*F2(i,j);
	}
      }
    }
  }
  fft_inverse(Kappa,kappa);
}

void computeVector(NumMatrix<Complex>& matrix, int phasefactor) {
  double phase;
  for (int i=0; i< matrix.getRows(); i++) {
    for (int j=0; j< matrix.getColumns(); j++) {
      phase = arg(matrix(i,j))/phasefactor;
      matrix(i,j) = abs(matrix(i,j))*(cos(phase) + I*sin(phase));
    }
  }
}

//compute gamma, F and G from kappa
// gamma will be compute from the potential via finite differences
// F and G are computed directly from kappa in Fourier space
// Solving gamma in Fourier space and solving G with finite differences creates artefacts!
void createMapsFromKappa(NumMatrix<double>& kappa, double dx, NumMatrix<Complex>& gamma, NumMatrix<Complex>& F, NumMatrix<Complex>& G) {
  int N = kappa.getRows();
  int J = kappa.getColumns();

  // since the solution of psi below is periodic along y axis
  // L is defined by the columns of kappa
  double L = dx*J;
  NumMatrix<double> psi(N,J);
  computePotentialFromConvergence(kappa,psi,dx,L);
  computeShear(psi,gamma,dx);

  NumMatrix<Complex> Kappa(N,J);
  fft_forward(kappa,Kappa);
  computeFlexion1Fourier(Kappa,F,dx);
  computeFlexion2Fourier(Kappa,G,dx);
}

