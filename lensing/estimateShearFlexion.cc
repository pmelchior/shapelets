#include <ShapeletsLensing.h>
#include <fstream.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

typedef complex<double> Complex;
const Complex I = Complex(0,1);

namespace ublas = boost::numeric::ublas;

/// ***** THIS DOES NOT WORK: BACON & GOLDBERG (2005) *******
// void estimateUnlensedCoeffs(NumMatrix<double>& coeffs, NumMatrix<double>& average, NumMatrix<double>& unlensed) {
//   int N = coeffs.getRows();
//   for (int i=0; i<N; i++) {
//     for (int j=0; j<N; j++) {
//       // consider only i+j = even coeffs
//       if ((i+j)%2 == 0)
// 	unlensed(i,j) = coeffs(i,j);
//       // the odd coeffs come from average
//       else
// 	unlensed(i,j) = average(i,j);
//     }
//   }
// }
// // the coeffs should be given as nxn matrices only!
// void estimateShearFlexion(NumMatrix<double>& average, NumMatrix<double>& coeffs, double beta, Complex& shear, Complex& F, Complex& G, double& chi2) {
//   // bring coeffs into shape of average
//   // average is always bigger or equally sized than coeffs
//   int deltaDim = average.getRows() - coeffs.getRows();
//   if (deltaDim > 0) coeffs.changeDimension(deltaDim);
//   // if average is smaller coeffs.getRows()+3, we have to expand them
//   // since the dimensions will be increased by 3 due to lensing below
//   if (deltaDim < 3) {
//     average.changeDimension(3-deltaDim);
//     coeffs.changeDimension(3-deltaDim);
//   }
//   int nmax = coeffs.getRows()-1;

//   // define estimator for unlensed image
//   NumMatrix<double> unlensed(nmax+1,nmax+1);
//   estimateUnlensedCoeffs(coeffs,average,unlensed);

//   // convert coeff matrices into coeff vectors
//   NumMatrix<int> nVector;
//   int nCoeffs = getNCoeffs(nmax);
//   makeNVector(nVector,nCoeffs,nmax);
//   NumVector<double> coeffVector(nCoeffs), averageVector(nCoeffs), unlensedVector(nCoeffs);
//   matrixMapping(coeffs,coeffVector,0,nVector,nCoeffs);
//   matrixMapping(average,averageVector,0,nVector,nCoeffs);
//   matrixMapping(unlensed,unlensedVector,0,nVector,nCoeffs);
//   // create Shapelet for manipulations
//   Point2D xcentroid(0,0);
//   ShapeletsImage *s = new ShapeletsImage(unlensed,beta,xcentroid);
  
//   // now construct sheared coeffs
//   NumMatrix<double> lensed;
//   Complex gamma;
//   // shear 1
//   gamma=1;
//   s->shear(gamma);
//   lensed  = s->getCartesianCoeffs();
//   deltaDim = nmax+1 - lensed.getRows();
//   lensed.changeDimension(deltaDim); 
//   NumVector<double> f1(nCoeffs);
//   matrixMapping(lensed,f1,0,nVector,nCoeffs);
//   f1 -= unlensedVector; // because shear does f' = (1+ gamma*S) f
//   // shear 2
//   gamma = I;
//   s->setCartesianCoeffs(unlensed);
//   s->shear(gamma);
//   lensed  = s->getCartesianCoeffs();
//   lensed.changeDimension(deltaDim); 
//   NumVector<double> f2(nCoeffs);
//   matrixMapping(lensed,f2,0,nVector,nCoeffs);
//   f2 -= unlensedVector;

//   // now the same for the flexion
//   NumMatrix<double> dGamma(2,2);
//   dGamma(0,0) = 1;
//   s->setCartesianCoeffs(unlensed);
//   s->flex(dGamma);
//   lensed  = s->getCartesianCoeffs();
//   deltaDim = nmax + 1 - lensed.getRows();
//   lensed.changeDimension(deltaDim);
//   NumVector<double> f3(nCoeffs);
//   matrixMapping(lensed,f3,0,nVector,nCoeffs);
//   f3 -= unlensedVector;
//   // flex 0,1
//   dGamma(0,0) = 0;
//   dGamma(0,1) = 1;
//   s->setCartesianCoeffs(unlensed);
//   s->flex(dGamma);
//   lensed  = s->getCartesianCoeffs();
//   lensed.changeDimension(deltaDim);
//   NumVector<double> f4(nCoeffs);
//   matrixMapping(lensed,f4,0,nVector,nCoeffs);
//   f4 -= unlensedVector;
//   // flex 1,0
//   dGamma(0,1) = 0;
//   dGamma(1,0) = 1;
//   s->setCartesianCoeffs(unlensed);
//   s->flex(dGamma);
//   lensed = s->getCartesianCoeffs();
//   lensed.changeDimension(deltaDim);
//   NumVector<double> f5(nCoeffs);
//   matrixMapping(lensed,f5,0,nVector,nCoeffs);
//   f5 -= unlensedVector;
//   // flex 1,1
//   dGamma(1,0) = 0;
//   dGamma(1,1) = 1;
//   s->setCartesianCoeffs(unlensed);
//   s->flex(dGamma);
//   lensed = s->getCartesianCoeffs();
//   lensed.changeDimension(deltaDim);
//   NumVector<double> f6(nCoeffs);
//   matrixMapping(lensed,f6,0,nVector,nCoeffs);
//   f6 -= unlensedVector;

//   // construct LS matrix
//   NumMatrix<double> Xt(6,nCoeffs);
//   for (int i=0; i< nCoeffs; i++) {
//     Xt(0,i) = -f1(i); // gamma1
//     Xt(1,i) = -f2(i); // gamma2
//     Xt(2,i) = -f3(i); // gamma1,1
//     Xt(3,i) = -f4(i); // gamma1,2
//     Xt(4,i) = -f5(i); // gamma2,1
//     Xt(5,i) = -f6(i); // gamma2,2
//   }
//   NumMatrix<double> LS = (Xt*Xt.transpose());
//   LS = LS.invert();
//   LS = LS*Xt;

//   // compute difference between lensed and unlensed coeffs
//   // this is mu - f in Goldberg/Bacon (2005) eq. 35
//   NumVector<double> delta = averageVector;
//   delta -= coeffVector;

//   // now compute the shear and flexions
//   NumVector<double> SF;
//   SF = LS * delta;
//   shear = SF(0) + I*SF(1);
//   F = SF(2) + SF(5) + I*(SF(4) - SF(3));
//   G = SF(2) - SF(5) + I*(SF(4) + SF(3));
  
//   NumVector<double> chi;
//   chi = (Xt.transpose())*SF;
//   chi += delta;
//   chi2 = chi*chi;
//   delete s;
// }

void computeAverageLSMatrixNVector(std::string averageFile, NumVector<double>& averageVector, NumMatrix<double>& LS, NumMatrix<double>& X, NumMatrix<int>& nVector, int& nmax) {
  NumMatrix<double> average;
  ShapeletObject *s = new ShapeletObject(averageFile);
  average = s->getCartesianCoeffs();
  
  // because of lensing coeff matrix will be expanded by 3
  // the averageVector has to reflect that
  average.resize(average.getRows()+3,average.getColumns()+3);
  nmax = average.getRows()-1;
  
  // convert coeff matrices into coeff vectors
  int nCoeffs = getNCoeffs(nmax);
  makeNVector(nVector,nCoeffs,nmax);
  averageVector = NumVector<double>(nCoeffs);
  matrixMapping(average,averageVector,0,nVector,nCoeffs);
  
  //// create Shapelet for manipulations
  //Point2D xcentroid(0,0);
  //ShapeletObject *s = new ShapeletObject(average,1,xcentroid);
  
  // now construct sheared coeffs
  NumMatrix<double> lensed;
  NumMatrix<double> Xt(6,nCoeffs); // the minimizing function matrix
  
  Complex gamma;
  // shear 1 -> first row of Xt matrix
  gamma=1;
  s->shear(gamma);
  lensed  = s->getCartesianCoeffs();
  //int deltaDim = nmax+1 - lensed.getRows();
  //lensed.changeDimension(deltaDim,deltaDim);
  lensed.resize(nmax+1,nmax+1);
  ublas::matrix_row< ublas::matrix<double> >mr (Xt,0);
  matrixMapping(lensed,mr,0,nVector,nCoeffs);
  mr -= averageVector;
  // shear 2
  gamma = I;
  s->setCartesianCoeffs(average);
  s->shear(gamma);
  lensed  = s->getCartesianCoeffs();
  //lensed.changeDimension(deltaDim,deltaDim); 
  lensed.resize(nmax+1,nmax+1);
  ublas::matrix_row< ublas::matrix<double> > mr1 (Xt,1);
  matrixMapping(lensed,mr1,0,nVector,nCoeffs);
  mr1 -= averageVector;

  // now the same for the flexion
  NumMatrix<double> dGamma(2,2);
  dGamma(0,0) = 1;
  s->setCartesianCoeffs(average);
  s->flex(dGamma);
  lensed  = s->getCartesianCoeffs();
  //deltaDim = nmax + 1 - lensed.getRows();
  //lensed.changeDimension(deltaDim,deltaDim);
  lensed.resize(nmax+1,nmax+1);
  ublas::matrix_row< ublas::matrix<double> >mr2 (Xt,2);
  matrixMapping(lensed,mr2,0,nVector,nCoeffs);
  mr2 -= averageVector;
  // flex 0,1
  dGamma(0,0) = 0;
  dGamma(0,1) = 1;
  s->setCartesianCoeffs(average);
  s->flex(dGamma);
  lensed  = s->getCartesianCoeffs();
  //lensed.changeDimension(deltaDim,deltaDim);
  lensed.resize(nmax+1,nmax+1);
  ublas::matrix_row< ublas::matrix<double> > mr3 (Xt,3);
  matrixMapping(lensed,mr3,0,nVector,nCoeffs);
  mr3 -= averageVector;
  // flex 1,0
  dGamma(0,1) = 0;
  dGamma(1,0) = 1;
  s->setCartesianCoeffs(average);
  s->flex(dGamma);
  lensed = s->getCartesianCoeffs();
  //lensed.changeDimension(deltaDim,deltaDim);
  lensed.resize(nmax+1,nmax+1);
  ublas::matrix_row< ublas::matrix<double> > mr4(Xt,4);
  matrixMapping(lensed,mr4,0,nVector,nCoeffs);
  mr4 -= averageVector;
  // flex 1,1
  dGamma(1,0) = 0;
  dGamma(1,1) = 1;
  s->setCartesianCoeffs(average);
  s->flex(dGamma);
  lensed = s->getCartesianCoeffs();
  //lensed.changeDimension(deltaDim,deltaDim);
  lensed.resize(nmax+1,nmax+1);
  ublas::matrix_row< ublas::matrix<double> > mr5 (Xt,5);
  matrixMapping(lensed,mr5,0,nVector,nCoeffs);
  mr5 -= averageVector;

  // construct LS matrix
  X = Xt.transpose();
  LS = (Xt*X);
  LS = LS.invert();
  LS = LS*Xt;

  // clean up
  delete s;
  //average.changeDimension(-3);
}

// the coeffs should be given as nxn matrices only!
void estimateShearFlexion(NumMatrix<double>& LS, NumMatrix<double>& X, const NumVector<double>& averageVector, const NumMatrix<int>& nVector, int nmax, NumMatrix<double>& coeffs, Complex& shear, Complex& F, Complex& G, double& chi2) {
  // bring coeffs into shape of average
  // average is always bigger or equally sized than coeffs
  
  //int deltaDim = nmax+1 - coeffs.getRows();
  //if (deltaDim > 0) coeffs.changeDimension(deltaDim,deltaDim);
  //coeffs.changeDimension(3,3);
  if (nmax+1 > coeffs.getRows())
    coeffs.resize(nmax+1,nmax+1);
  coeffs.resize(coeffs.getRows()+3, coeffs.getColumns()+3);

  int nCoeffs = getNCoeffs(nmax);
  NumVector<double> coeffVector(nCoeffs);
  matrixMapping(coeffs,coeffVector,(bool)0,nVector,nCoeffs);
  
  // compute difference between lensed and unlensed coeffs
  // this is mu - f in Goldberg/Bacon (2005) eq. 35
  NumVector<double> delta = coeffVector;
  delta -= averageVector;

  // now compute the shear and flexions
  NumVector<double> SF;
  SF = LS * delta;
  shear = SF(0) + I*SF(1);
  F = SF(2) + SF(5) + I*(SF(4) - SF(3));
  G = SF(2) - SF(5) + I*(SF(4) + SF(3));
  
  NumVector<double> chi;
  chi = X*SF;
  chi -= delta;
  chi2 =  chi*chi;

  coeffs.resize(coeffs.getRows()-3, coeffs.getColumns()-3);
}



// calculate the mean and the SD of the mean
// for gamma, F and G with weights coming from the minimization chi2
void getMeanSigmaMean(NumVector<Complex>& gammas, NumVector<Complex>& Fs, NumVector<Complex>& Gs, NumVector<double>& chi2, Complex& gamma, Complex& F, Complex& G, Complex& sigma_gamma, Complex& sigma_F, Complex& sigma_G, int NOBJ) {
  NumVector<Complex> mu_gamma(NOBJ), mu_F(NOBJ), mu_G(NOBJ), weight(NOBJ);
  // this computes the series of means of the sample at length n
  for (int n=0; n < NOBJ; n++) {
    mu_gamma(n) = gammas(n)/chi2(n);
    mu_F(n) = Fs(n)/chi2(n);
    mu_G(n) = Gs(n)/chi2(n);
    weight(n) = 1./chi2(n);
    if (n>0) {
      mu_gamma(n) += mu_gamma(n-1)*weight(n-1);
      mu_F(n) += mu_F(n-1)*weight(n-1);
      mu_G(n) += mu_G(n-1)*weight(n-1);
      weight(n) += weight(n-1);
    }
    mu_gamma(n) /= weight(n);
    mu_F(n) /= weight(n);
    mu_G(n) /= weight(n);
  }
  gamma = mu_gamma(NOBJ-1);
  F = mu_F(NOBJ-1);
  G = mu_G(NOBJ-1);
  
    // calculate the the standard deviation of the means
  sigma_gamma = sigma_F = sigma_G = 0;
  for (int n=0; n < NOBJ; n++) {
    sigma_gamma += gsl_pow_2(real(mu_gamma(n)) - real(gamma));
    sigma_gamma += I*gsl_pow_2(imag(mu_gamma(n)) - imag(gamma));
    sigma_F += gsl_pow_2(real(mu_F(n)) - real(F));
    sigma_F += I*gsl_pow_2(imag(mu_F(n)) - imag(F));
    sigma_G += gsl_pow_2(real(mu_G(n)) - real(G));
    sigma_G += I*gsl_pow_2(imag(mu_G(n)) - imag(G));
  }
  sigma_gamma = sqrt(real(sigma_gamma)/(NOBJ-1)) + I*sqrt(imag(sigma_gamma)/(NOBJ-1));
  sigma_F = sqrt(real(sigma_F)/(NOBJ-1)) + I*sqrt(imag(sigma_F)/(NOBJ-1));
  sigma_G = sqrt(real(sigma_G)/(NOBJ-1)) + I*sqrt(imag(sigma_G)/(NOBJ-1));
}

// computes the average and sigma values of shear and flexions 
// for all objects inside one image pixel
// The sigma0 are the errors in the measurement coming from the scatter in 
// galaxy shapes
void estimateLensingInsidePixel(std::string path, int pixelX, int pixelY, int NOBJ, NumMatrix<double>& LS, NumMatrix<double>& X, NumVector<double>& averageVector,NumMatrix<int>& nVector, int nmax, Complex& gamma, Complex& F, Complex& G, Complex& sigma_gamma, Complex& sigma_F, Complex& sigma_G) {
  NumVector<Complex> gammas(NOBJ), Fs(NOBJ), Gs(NOBJ);
  NumVector<double> chi2(NOBJ);
  for (int n=0; n < NOBJ; n++) {
    Complex gamma_, F_, G_;
    std::ostringstream filename;
    filename << path << "lensed_"<<pixelX<<"_"<<pixelY<<"_"<<n<<".sif";
    ShapeletObject *s = new ShapeletObject(filename.str());
    NumMatrix<double> coeffs = s->getCartesianCoeffs();
    estimateShearFlexion(LS,X,averageVector,nVector,nmax,coeffs,gammas(n),Fs(n),Gs(n),chi2(n));
    delete s;
  }
  getMeanSigmaMean(gammas,Fs,Gs,chi2,gamma,F,G,sigma_gamma,sigma_F,sigma_G,NOBJ);

  // this computes the prediction of the lensing parameters
  // in dependence of the number of galaxies per pixel
//   int n=1;
//   while (n < NOBJ) {
//     getMeanSigmaMean(gammas,Fs,Gs,chi2,gamma,F,G,sigma_gamma,sigma_F,sigma_G,n);
//     std::cout << n << " " << real(gamma) << " " << real(sigma_gamma) << " ";
//     std::cout << imag(gamma) << " " << imag(sigma_gamma) << " ";
//     std::cout << real(F) << " " << real(sigma_F) << " ";
//     std::cout << imag(F) << " " << imag(sigma_F) << " ";
//     std::cout << real(G) << " " << real(sigma_G) << " ";
//     std::cout << imag(G) << " " << imag(sigma_G) << std::endl;
//     if (n < 50) n+=2;
//     else if (n < 200) n+=10;
//     else n+=100;
//   }
}

// compute error on shear and flexion measurement coming from the scatter
// in galaxy shapes that are compared to the average galaxy shape
// the mean here should be 0 within the errors as long the galaxies were unlensed.
void estimateInstrinsicLensing(std::string path, std::string averageFile, int NOBJ, Complex& gamma, Complex& F, Complex& G, Complex& sigma_gamma, Complex& sigma_F, Complex& sigma_G) {
  // compute the vector representation of the averageCoeffs
  // and the LS matrix based on those coeffs
  NumVector<double> averageVector;
  NumMatrix<double> LS,X;
  NumMatrix<int> nVector;
  int nmax;
  computeAverageLSMatrixNVector(averageFile,averageVector,LS,X,nVector,nmax);

  NumVector<Complex> gammas(NOBJ), Fs(NOBJ), Gs(NOBJ);
  NumVector<double> chi2(NOBJ);
  for (int n=0; n < NOBJ; n++) {
    Complex gamma_, F_, G_;
    std::ostringstream filename;
    filename << path << "shapelets_" << n << ".sif";
    ShapeletObject *s = new ShapeletObject(filename.str());
    NumMatrix<double> coeffs = s->getCartesianCoeffs();
    estimateShearFlexion(LS,X,averageVector,nVector,nmax,coeffs,gammas(n),Fs(n),Gs(n),chi2(n));
    delete s;
  }

//   gsl_histogram2d * h = gsl_histogram2d_alloc (20,20);
//   gsl_histogram2d_set_ranges_uniform (h,-0.05, 0.05, -0.05, 0.05);
//   for (int n=0; n < NOBJ; n++)
//     gsl_histogram2d_increment (h, real(Gs(n)),imag(Gs(n)));
//   gsl_histogram2d_fprintf (stdout, h, "%g", "%g");
//   gsl_histogram2d_free (h);
  
  getMeanSigmaMean(gammas,Fs,Gs,chi2,gamma,F,G,sigma_gamma,sigma_F,sigma_G,NOBJ);
}

void estimateLensingMaps(std::string path, std::string averageFile, int NOBJ, NumMatrix<Complex>& gamma, NumMatrix<Complex>& F, NumMatrix<Complex>& G, NumMatrix<Complex>& sigma_gamma, NumMatrix<Complex>& sigma_F, NumMatrix<Complex>& sigma_G) {
  // now compute the vector representation of the averageCoeffs
  // and the LS matrix based on those coeffs
  // This has to be done only once!
  NumVector<double> averageVector;
  NumMatrix<double> LS,X;
  NumMatrix<int> nVector;
  int nmax;
  computeAverageLSMatrixNVector(averageFile,averageVector,LS,X,nVector,nmax);
  // now run over each pixel in th image
  for (int i=0; i < gamma.getRows(); i++) {
    for (int j=0; j< gamma.getColumns(); j++) { 
      estimateLensingInsidePixel(path,i,j,NOBJ,LS,X,averageVector,nVector,nmax,gamma(i,j),F(i,j),G(i,j),sigma_gamma(i,j), sigma_F(i,j), sigma_G(i,j));
    }
  }
}
