#ifndef SHAPELETS_LENSING_ANALYSIS_H
#define SHAPELETS_LENSING_ANALYSIS_H
#include <ShapeletObject.h>

// FIXME: set up makefile, create one central source file with the pipeline.


/// functions from FFT.cc
void fft_forward(NumVector<double>& f, NumVector<complex<double> >& F);
void fft_forward(NumVector<complex<double> >& f, NumVector<complex<double> >& F);
void fft_forward(NumMatrix<double>& f, NumMatrix<complex<double> >& F);
void fft_forward(NumMatrix<complex<double> >& f, NumMatrix<complex<double> >& F);
void fft_sin_forward(NumVector<double>& f, NumVector<double>& F);
void fft_cos_forward(NumVector<double>& f, NumVector<double>& F);
void fft_inverse(NumVector<complex<double> >& F, NumVector<double>& f);
void fft_inverse(NumVector<complex<double> >& F, NumVector<complex<double> >& f);
void fft_inverse(NumMatrix<complex<double> >& F, NumMatrix<double>& f);
void fft_inverse(NumMatrix<complex<double> >& F, NumMatrix<complex<double> >& f);
void fft_sin_inverse(NumVector<double>& F, NumVector<double>& f);
void fft_cos_inverse(NumVector<double>& F, NumVector<double>& f);
double wavenumber(int j, double dx, int N);

/// the functions of Maps.cc
void computePotentialFromConvergence(NumMatrix<double>& kappa, NumMatrix<double>& psi, double dx, double L);
void computeShearFourier(NumMatrix<complex<double> >& Kappa, NumMatrix<complex<double> >& shear, double dx);
void computeFlexion1Fourier(NumMatrix<complex<double> >& Kappa, NumMatrix<complex<double> >& F, double dx);
void computeFlexion2Fourier(NumMatrix<complex<double> >& Kappa, NumMatrix<complex<double> >& G, double dx);
void computeConvergence(NumMatrix<double>& psi, NumMatrix<complex<double> >& kappa, double dx);
void computeShear(NumMatrix<double>& psi, NumMatrix<complex<double> >& shear, double dx);
void computeFlexion(NumMatrix<double>& psi, bool number, NumMatrix<complex<double> >& flexion, double dx);
void computeConvergenceFromShear(NumMatrix<complex<double> >& shear, NumMatrix<complex<double> >& kappa, double L, double dx);
void computeConvergenceFromFlexion1(NumMatrix<complex<double> >& flexion,NumMatrix<complex<double> >& kappa, double L, double dx);
void computeConvergenceFromFlexionFourier(NumMatrix<complex<double> >& flexion, bool number, NumMatrix<complex<double> >& kappa, double dx);
void computeVector(NumMatrix<complex<double> >& matrix, int phasefactor);
void createMapsFromKappa(NumMatrix<double>& kappa, double dx, NumMatrix<complex<double> >& gamma, NumMatrix<complex<double> >& F, NumMatrix<complex<double> >& G); 

/// functions from createShapeletsImage.cc
/// the images created here are all flux-normalized.
void createShapeletImages(NumMatrix<double>& averageCoeffs, NumMatrix<double>& sigmaCoeffs, double betamin, double betamax, std::string path);
void createShapeletImagesPCA( double betamin, double betamax, std::string path, int N);
void averageShapeletCoeffs(NumMatrix<double>& average, double& beta, std::string listfile, bool normalize);
void createLensedShapeletImages(double dx, NumMatrix<double>& kappa, NumMatrix<complex<double> >& gamma, NumMatrix<complex<double> >& F, NumMatrix<complex<double> >& G, std::string listfilename, std::string writeDirectory, int NOBJ);


/// functions from estimateShearFlexion.cc
/// is is crucial, that the lensed galaxies have been flux-normalized and set to beta = 1
/// otherwise the minimization is dominated by noise due to missing comparability
/// - flux difference would introduce factor between mu and f
/// - beta difference would change lensing behaviour for flexion since it depends on 
///   second brightness moments
void estimateInstrinsicLensing(std::string path, std::string averageFile, int NOBJ,complex<double>& gamma, complex<double>& F, complex<double>& G, complex<double>& sigma_gamma, complex<double>& sigma_F, complex<double>& sigma_G);
void estimateLensingInsidePixel(std::string path, int pixelX, int pixelY, int NOBJ, NumMatrix<double>& LS, NumMatrix<double>& X, NumVector<double>& averageVector,NumMatrix<int>& nVector, int nmax, complex<double>& gamma, complex<double>& F, complex<double>& G, complex<double>& sigma_gamma, complex<double>& sigma_F, complex<double>& sigma_G);
void estimateShearFlexion(NumMatrix<double>& LS, NumMatrix<double>& X, NumVector<double>& averageVector,NumMatrix<int>& nVector, int nmax, NumMatrix<double>& coeffs, complex<double>& shear, complex<double>& F, complex<double>& G, double& chi2);
void estimateLensingMaps(std::string path, std::string averageFile, int NOBJ, NumMatrix<complex<double> >& gamma, NumMatrix<complex<double> >& F, NumMatrix<complex<double> >& G, NumMatrix<complex<double> >& sigma_gamma, NumMatrix<complex<double> >& sigma_F, NumMatrix<complex<double> >& sigma_G);
void computeAverageLSMatrixNVector(std::string averageFile, NumVector<double>& averageVector, NumMatrix<double>& LS, NumMatrix<double>& X, NumMatrix<int>& nVector, int& nmax);

#endif
