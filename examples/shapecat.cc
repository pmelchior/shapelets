#include <ShapeLens.h>
#include <tclap/CmdLine.h>
#include <fstream>

void ellipticity(const NumMatrix<data_t>& Q, data_t& e1, data_t& e2, data_t& e, data_t& theta) {
  complex<data_t> I(0,1);
  complex<data_t> Q11(Q(0,0),0),Q22(Q(1,1),0),Q12(Q(0,1),0);
  complex<data_t> denom = Q11+Q22;// + 2.*sqrt(Q11*Q22-Q12*Q12);
  complex<data_t> epsilon = (Q11 - Q22 + data_t(2)*I*Q12)/denom;
  
  e1 = real(epsilon);
  e2 = imag(epsilon);
  e = sqrt(e1*e1 + e2*e2);
  data_t r = (1-e)/(1+e);
  if (Q(0,0)-Q(1,1)!=0)
    theta = 180./M_PI*(0.5*atan(2*Q(1,0)/(Q(0,0)-Q(1,1))));
  else
    theta = 0;
  // this ensures that theta has the same sign as Q12
  if ((theta < 0 && Q(0,1) > 0) || (theta > 0 && Q(0,1) < 0))
    theta *= -1;
}

int main(int argc, char *argv[]) {
  TCLAP::CmdLine cmd("Compute various shape statistic based on shapelet coefficients", ' ', "0.3");
  TCLAP::ValueArg<std::string> listArg("l","list","File that lists the SIF files",true,"","string", cmd);
  TCLAP::ValueArg<std::string> outputArg("o","output","Name of output file for shape catalog",false,"","string", cmd);
  cmd.parse( argc, argv );

  // open ensembles of SIF files
  ShapeletObjectList sl(listArg.getValue());

  std::ofstream output;
  if (outputArg.isSet()) {
    output.open(outputArg.getValue().c_str(),std::ios::out);
    output << "# ID BETA NMAX CHI2 R_S FLUX XC1 XC2 Q11 Q12 Q22 E1 E2 ELLIP THETA RMS FLAGS" << std::endl;
  } else
    std::cout << "# ID BETA NMAX CHI2 R_S FLUX XC1 XC2 Q11 Q12 Q22 E1 E2 ELLIP THETA RMS FLAGS" << std::endl;
  
  int nmax;
  data_t beta, chi2, flux, e1,e2,e,theta, RMS, Rs;
  Point2D scentroid;
  NumMatrix<data_t> Q(2,2);
  for (ShapeletObjectList::iterator iter = sl.begin(); iter != sl.end() ; iter++) {
    // get shapelet parameters used in decomposition
    beta = (*iter)->getBeta();
    nmax = (*iter)->getNMax();
    chi2 = (*iter)->getChiSquare();
    Rs = sqrt(((*iter)->getCoeffs())*((*iter)->getCoeffs()));

    // get flux and centroid, 2nd moments and RMS radius
    flux = (*iter)->getShapeletFlux();
    scentroid = (*iter)->getShapeletCentroid();
    Q = (*iter)->getShapelet2ndMoments();
    RMS = (*iter)->getShapeletRMSRadius();
    // compute ellipticity and orientation from Q
    ellipticity(Q,e1,e2,e,theta);

    // output statistics
    if (outputArg.isSet())
      output << (*iter)->getObjectID() << " " << beta << " " << nmax << " " << chi2 << " " << Rs << " " << flux << " " << scentroid(0) << " " << scentroid(1) << " " << Q(0,0) << " " << Q(0,1) << " " << Q(1,1) << " " << e1 << " " << e2 << " " << e << " " << theta << " " << RMS  << " " << (*iter)->getFlags().to_string<char,std::char_traits<char>,std::allocator<char> >() << std::endl;
    else
      std::cout << (*iter)->getObjectID() << " " << beta << " " << nmax << " " << chi2 << " " <<  Rs << " " << flux << " " << scentroid(0) << " " << scentroid(1) << " " << Q(0,0) << " " << Q(0,1) << " " << Q(1,1) << " " << e1 << " " << e2 << " " << e << " " << theta << " " << RMS << " " << (*iter)->getFlags().to_string<char,std::char_traits<char>,std::allocator<char> >() << std::endl;
    
  }
  if (outputArg.isSet())
    output.close();
}
  
