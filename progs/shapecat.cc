#include <shapelens/ShapeLens.h>
#include <tclap/CmdLine.h>
#include <fstream>

using namespace shapelens;

void ellipticity(Quadrupole& Q, data_t& e1, data_t& e2, data_t& e, data_t& theta) {
  complex<data_t> I(0,1);
  complex<data_t> Q11(Q(0,0),0),Q22(Q(1,1),0),Q12(Q(0,1),0);
  complex<data_t> denom = Q11+Q22 + 2.*sqrt(Q11*Q22-Q12*Q12);
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
  TCLAP::CmdLine cmd("Compute various shape statistic based on shapelet coefficients", ' ', "0.5");
  TCLAP::ValueArg<std::string> listArg("l","list","File that lists the SIF files",true,"","string");
  TCLAP::ValueArg<std::string> input("i","input","SIF file to analyze",true,"","string");
  TCLAP::ValueArg<std::string> table("t","table","Table in DB to load SObjs (uses where is present)",true,"","string");
  std::vector<TCLAP::Arg*> inputs;
  inputs.push_back(&listArg);
  inputs.push_back(&input);
  inputs.push_back(&table);
  cmd.xorAdd(inputs);
  TCLAP::ValueArg<std::string> where("w","where","Where statement to select SObjs in table",false,"","string",cmd);
  TCLAP::ValueArg<std::string> kernel("k","kernel","Deconvolve from this kernel SIF file",false,"","string",cmd);
  TCLAP::ValueArg<unsigned int> truncate("T","truncate","Truncate coefficient set at given order n_max",false,0,"unsigned int",cmd);
  TCLAP::ValueArg<unsigned int> diamond("D","diamond","Truncate coefficient set to diamond shape of given order n_max",false,0,"unsigned int",cmd);
  TCLAP::ValueArg<std::string> outputArg("o","output","Name of output file for shape catalog",false,"","string", cmd);
  cmd.parse( argc, argv );

  // open ensembles of SIF files
  ShapeletObjectList sl;
  if (listArg.isSet())
    sl = ShapeletObjectList(listArg.getValue());
  else if (input.isSet())
    sl.push_back(boost::shared_ptr<ShapeletObject>(new ShapeletObject(input.getValue())));
  else {
    ShapeletObjectDB db;
    db.selectTable(table.getValue());
    sl = db.load(where.getValue());
  }

  ShapeletObject sk;
  if (kernel.isSet())
    sk = ShapeletObject(kernel.getValue());
  
  std::ofstream output;
  if (outputArg.isSet()) {
    output.open(outputArg.getValue().c_str(),std::ios::out);
    output << "# ID BETA NMAX CHI2 R_S FLUX XC1 XC2 Q11 Q12 Q22 E1 E2 ELLIP THETA RMS FLAGS" << std::endl;
  } else
    std::cout << "# ID BETA NMAX CHI2 R_S FLUX XC1 XC2 Q11 Q12 Q22 E1 E2 ELLIP THETA RMS FLAGS" << std::endl;
  
  int nmax;
  data_t beta, chi2, flux, e1,e2,e,theta, RMS, Rs;
  Point2D<data_t> scentroid;
  Quadrupole Q;
  for (ShapeletObjectList::iterator iter = sl.begin(); iter != sl.end() ; iter++) {
    // deconvolve if demanded
    if (kernel.isSet())
      (*iter)->deconvolve(sk.getCoeffs(), sk.getBeta());

    // apply truncation of demanded
    if(truncate.isSet()) {
      (*iter)->setNMax(truncate.getValue());
    }
    else if (diamond.isSet()) {
      CoefficientVector<complex<data_t> > pcoeffs = (*iter)->getPolarCoeffs();
      int limit;
      if (diamond.getValue() % 2 == 0)
	limit = diamond.getValue();
      else
	limit = diamond.getValue() + 1;
      // lower order
      pcoeffs.setNMax(diamond.getValue());
      for (int n = limit/2; n <= diamond.getValue(); n++)
	for (int m = -n; m <= n; m+=2)
	  if (abs(m) > limit-n)
	    pcoeffs(n,m) = complex<data_t>(0,0);
      (*iter)->setPolarCoeffs(pcoeffs);
    }

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
  
