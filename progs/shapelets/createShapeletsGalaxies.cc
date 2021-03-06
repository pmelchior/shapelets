#include <shapelens/ShapeLens.h>
#include <shapelens/shapelets/ShapeletObject.h>
#include <shapelens/shapelets/ShapeletObjectList.h>
#include <fstream>
#include <string>
#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <list>
#include <tclap/CmdLine.h>

using namespace shapelens;
typedef std::complex<data_t> Complex;

/// Creates \f$N\f$ galaxy images in shapelet space.
/// The output is stored as "shapelets_<i>.sif", where \f$0 \le i < N\f$
/// in the directory specified by path.
/// The coeffs are build by gaussian distributions with width sigma arround
/// the average.
/// A file "shapelets.ls" and a sif file "average.sif" with the average coeffs 
/// and beta derived from this sample of galaxies will be stored also in 
/// the directory.\n
/// Note: Don't forget the / at the end of the path.
void createShapeletImages(CoefficientVector<data_t>& averageCoeffs, CoefficientVector<data_t>& sigmaCoeffs, data_t betamin, data_t betamax, std::string path, int N) {
  // create directory if it doesn't exist
  std::ostringstream dircall;
  dircall << "mkdir -p " << path << "&> /dev/null";
  std::system(dircall.str().c_str());

  //create a generator chosen by the 
  // environment variable GSL_RNG_TYPE
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  // build template shapelet
  Point<data_t> xcentroid(0,0);
  CoefficientVector<data_t> coeffs(averageCoeffs.getNMax());
  data_t beta = 1;
  ShapeletObject s (coeffs,beta,xcentroid);
  
   // now different galaxies, all normalized and with beta = 1
  std::ostringstream sif_filename;
  std::ofstream listfile((path + "shapelets.ls").c_str());
  for (int n = 0; n < N; n++) {
    // add gaussian scatter on coeffs
    for (int i=0; i<coeffs.size(); i++)
      coeffs(i) = averageCoeffs(i) + gsl_ran_gaussian (r, sigmaCoeffs(i));

    s.setCoeffs(coeffs);

    // beta: uniform between betamin and betamax
    beta =  betamin + gsl_rng_uniform(r)*(betamax-betamin);
    s.setBeta(beta);
    // set grid
    int range = std::max((int)ceil(15*beta),40);
    if (range%2==1) range++;
    Grid grid = Grid(-range/2,-range/2,range,range);
    s.setGrid(grid);
        
    // since images are aligned along x axis
    // rotate them by any angle between 0 and 2 pi
    s.rotate(2*M_PI*gsl_rng_uniform(r));
    //s.shear(std::complex<data_t>(0.2,0));
    // normalize flux for comparability later
    s.brighten(1./s.getShapeletFlux());
    s.brighten(beta*beta);

    // save as sif file
    sif_filename.str("");
    sif_filename << path << "shapelets_" << n << ".sif";
    listfile << sif_filename.str() << std::endl;
    s.setHistory("");
    s.save(sif_filename.str());
  }
  gsl_rng_free (r);
  listfile.close();

  // create average sif file
  CoefficientVector<data_t> mean, std_mean;
  data_t averageBeta;
  std::ostringstream filename;
  filename << path << "shapelets.ls";
  ShapeletObjectList sl(filename.str());
  sl.average(mean, std_mean, averageBeta);
  ShapeletObject a(mean,averageBeta,xcentroid);
  a.setErrors(std_mean);
  
  filename.str("");
  filename << path << "average.sif";
  a.save(filename.str());
}

/// Creates shapelet image base on the dominant shapelet states
/// and their scatter, derived in Kelly & McKay (2003): astro-ph/0307395.
/// \f$N\f$ denotes the number of created galaxies.
void createShapeletImagesPCA(data_t betamin, data_t betamax,std::string path, int N) {
  NumMatrix<data_t> coeffs(9,9), sigma(9,9);
  // this is the table 1 in the reference
  coeffs(0,0) = 0.293;
  coeffs(4,0) = 0.047;
  coeffs(2,0) = 0.059;
  coeffs(0,4) = 0.026;
  coeffs(6,0) = 0.022;
  coeffs(2,2) = 0.021;
  coeffs(0,2) = 0.030;
  coeffs(8,0) = 0.018;
  coeffs(0,1) = 0.016;
  coeffs(1,0) = 0.014;
  // this comes also from table 1
  sigma(0,0) = 0.076;
  sigma(4,0) = 0.012;
  sigma(2,0) = 0.021;
  sigma(0,4) = 0.011;
  sigma(6,0) = 0.010;
  sigma(2,2) = 0.010;
  sigma(0,2) = 0.020;
  sigma(8,0) = 0.009;
  sigma(0,1) = 0.015;
  sigma(1,0) = 0.012;
  // now adding additional scatter in minor coeffs
  sigma(0,8) = 0.009; // copy from sigma(8,0)
  sigma(0,6) = 0.010; // copy from sigma(6,0)
  sigma(1,1) = 0.003;
  sigma(3,3) = 0.005;
  sigma(4,4) = 0.006;
  sigma(0,3) = sigma(3,0) = 0.005;
  sigma(0,5) = sigma(5,0) = 0.005;
  sigma(0,7) = sigma(7,0) = 0.005;
  sigma(1,2) = sigma(2,1) = 0.003;
  sigma(1,3) = sigma(3,1) = 0.003;
  sigma(1,4) = sigma(4,1) = 0.003;
  sigma(1,5) = sigma(5,1) = 0.003;
  sigma(1,6) = sigma(6,1) = 0.003;
  sigma(1,7) = sigma(7,1) = 0.003;
  sigma(2,6) = sigma(6,2) = 0.002;
  sigma(2,3) = sigma(3,2) = 0.002;
  sigma(2,4) = sigma(4,2) = 0.002;
  sigma(2,5) = sigma(5,2) = 0.002;
  sigma(3,4) = sigma(4,3) = 0.001;
  sigma(3,5) = sigma(5,3) = 0.001;
  CoefficientVector<data_t> mean(coeffs), std(sigma);
  createShapeletImages(mean,std,betamin,betamax,path,N);
}

void createLensedShapeletImages(data_t dx, NumMatrix<data_t>& kappa, NumMatrix<Complex>& gamma, NumMatrix<Complex>& F, NumMatrix<Complex>& G, std::string listfilename, std::string writeDirectory, int NOBJ) {
  int N = kappa.getRows();
  int J = kappa.getColumns(); 
  std::string filename;

  // create directory if it doesn't exist
  std::ostringstream dircall;
  dircall << "mkdir -p " << writeDirectory << "&> /dev/null";
  std::system(dircall.str().c_str());

  // write applied lensing parameters to file
  std::ofstream mapFile;
  std::string mapFilename;
  mapFilename = writeDirectory + "map.dat";
  mapFile.open(mapFilename.c_str());
  for (int i=0; i< N; i++) {
    for (int j=0; j< J; j++) {
      mapFile << i*dx << " " << j*dx << " " << kappa(i,j) << " ";
      mapFile << real(gamma(i,j)) << " " << imag(gamma(i,j)) << " ";
      mapFile << real(F(i,j)) << " " << imag(F(i,j)) << " ";
      mapFile << real(G(i,j)) << " " << imag(G(i,j)) << " ";
      mapFile << std::endl;
    }
  }   
  mapFile.close();

  // read in filelist of unlensed sif file object
  std::vector<std::string> filelist;
  std::ifstream fileWithList;
  fileWithList.open(listfilename.c_str());
  while (1) {
    fileWithList >> filename;
    filelist.push_back(filename);
    if (!fileWithList.good()) break;
  }
  fileWithList.close();
  int filenumber = filelist.size(), randomnumber;

  // open file for writing selection unlensed -> lensed to file selection.dat
  std::ostringstream selectionfilename;
  selectionfilename << writeDirectory << "selection.dat";
  std::ofstream selectionFile;
  selectionFile.open(selectionfilename.str().c_str());

  // for each point in the map select N_obj objects from sif files randomly
  // and apply lensing operations on them
  // store them as lensed_i_j_n.sif
  NumMatrix<data_t> dGamma(2,2);
  for (int i = 0; i < N; i++) {
    for (int j =0 ; j < J; j++) {
      for (int n = 0; n < NOBJ; n++) {
	// random number from filelist
	randomnumber = (int)floor((data_t)filenumber*rand()/(RAND_MAX+1.0));
	filename = filelist[randomnumber];
	ShapeletObject* s = new ShapeletObject(filename);
	// since mass sheat degeneracy: don't include kappa
	s->lens(0,gamma(i,j),F(i,j),G(i,j));
	std::ostringstream lensedfilename, selectionMessage;
	lensedfilename << writeDirectory << "lensed_" << i <<"_"<<j<<"_"<<n<<".sif";
	selectionMessage << n + j*NOBJ + i*J*NOBJ << ": " << filename << " -> " << lensedfilename.str() << std::endl;
	s->setHistory("");
	s->save(lensedfilename.str());
	selectionFile << selectionMessage.str();
	delete s;
      }
    }
  }
  selectionFile.close();
} 


int main(int argc, char *argv[]) {
  TCLAP::CmdLine cmd("Create set of simulated galaxy SIF files", ' ', "0.1");
  TCLAP::ValueArg<data_t> minArg("m","min","min(beta)",true,1,"data_t", cmd);
  TCLAP::ValueArg<data_t> maxArg("M","max","max(beta)",true,10,"data_t", cmd);
  TCLAP::ValueArg<int> NArg("N","NOBJ","total number of galaxies",true,100,"int", cmd);
  TCLAP::ValueArg<std::string> pathArg("p","path","path to store SIF files",true,"","string", cmd);
  cmd.parse(argc,argv);

  createShapeletImagesPCA(minArg.getValue(),maxArg.getValue(),pathArg.getValue(), NArg.getValue());
}
