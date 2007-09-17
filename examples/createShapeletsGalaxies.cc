#include <ShapeLens.h>
#include <fstream.h>
#include <string>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <list>

typedef complex<double> Complex;

/// Computes the average coeffs and the standard deviation of the average
// and beta from sif files listed in listfile
void averageShapeletCoeffs(NumMatrix<double>& average, NumMatrix<double>& std_mean, double& beta, std::string listfile) {
  //create new average and std_mean
  average.clear();
  std_mean.clear();

  // set up list of coeff matrices
  std::list<NumMatrix<double> > matrixList;
  beta = 0;
  int count = 0;
  std::ifstream files;
  std::string filename;
  files.open(listfile.c_str());
  while (1) {
    files >> filename;
     if (!files.good()) break;
    
    ShapeletObject* s = new ShapeletObject(filename);
    const NumMatrix<double>& coeffs = s->getCartesianCoeffs();
    matrixList.push_back(coeffs);
    // if new coeff matrix is bigger than current average matrix
    // expand average
    if (average.getRows() < coeffs.getRows() || average.getColumns() < coeffs.getColumns()) {
      average.resize_clear(coeffs.getRows(),coeffs.getColumns());
      std_mean.resize_clear(coeffs.getRows(),coeffs.getColumns());
    }
    beta += s->getBeta();
    count++;
    delete s;
  }
  files.close();
  
  // compute average beta
  beta /= count; // that's the number of considered images
  
  // now average over all coeffs and all matrices
  NumVector<double> entries(count);
  std::list<NumMatrix<double> >::iterator iter;
  for (int i=0; i<average.getRows(); i++) {
    for (int j=0; j<average.getColumns(); j++) {
      entries.clear();
      int n=0;
      for(iter = matrixList.begin(); iter != matrixList.end(); iter++ ) {
	if ((*iter).getRows() > i && (*iter).getColumns() > j) 
	  entries(n) = (*iter)(i,j);
	n++;
      }
      average(i,j) = entries.mean();
      std_mean(i,j) = entries.std()/sqrt(1.*count);
    }
  }
}

/// Creates \f$N\f$ galaxy images in shapelet space.
/// The output is stored as "shapelets_<i>.sif", where \f$0 \le i < N\f$
/// in the directory specified by path.
/// The coeffs are build by gaussian distributions with width sigma arround
/// the average.
/// A file "shapelets.ls" and a sif file "average.sif" with the average coeffs 
/// and beta derived from this sample of galaxies will be stored also in 
/// the directory.\n
/// Note: Don't forget the / at the end of the path.
void createShapeletImages(NumMatrix<double>& averageCoeffs, NumMatrix<double>& sigmaCoeffs, double betamin, double betamax, std::string path, int N) {
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
  Point2D xcentroid(0,0);
  NumMatrix<double> coeffs(averageCoeffs.getRows(),averageCoeffs.getColumns());
  double beta = 1;
  ShapeletObject *s = new ShapeletObject(coeffs,beta,xcentroid);
  
   // now different galaxies, all normalized and with beta = 1
  for (int n = 0; n < N; n++) {
    // add gaussian scatter on coeffs
    for (int i=0; i<coeffs.getRows();i++)
      for (int j=0; j<coeffs.getColumns();j++)
	coeffs(i,j) = averageCoeffs(i,j) + gsl_ran_gaussian (r, sigmaCoeffs(i,j));

    s->setCartesianCoeffs(coeffs);

    // beta: uniform between betamin and betamax
    beta =  betamin + gsl_rng_uniform(r)*(betamax-betamin);
    s->setBeta(beta);
    // set grid
    int range = GSL_MAX_INT((int)ceil(15*beta),40);
    if (range%2==1) range++;
    Grid grid = Grid(-range/2,range/2,1,-range/2,range/2,1);
    s->setGrid(grid);
        
    // since images are aligned along x axis
    // rotate them by any angle between 0 and 2 pi
    s->rotate(2*M_PI*gsl_rng_uniform(r));

    // normalize flux for comparability later
    s->brighten(1./s->integrate());
    s->brighten(beta*beta);

    // save as sif file
    std::ostringstream sif_filename;
    sif_filename << path << "shapelets_" << n << ".sif";
    s->setHistory("");
    s->save(sif_filename.str());
  }
  delete s;
  gsl_rng_free (r);

  // create list of files
  std::ostringstream job;
  job << "find " << path << " -type f -name 'shapelets_*.sif' > " << path << "shapelets.ls";
  std::system(job.str().c_str());

  // create average sif file
  NumMatrix<double> average, std_mean;
  double averageBeta;
  std::ostringstream filename;
  filename << path << "shapelets.ls";
  averageShapeletCoeffs(average,std_mean,averageBeta,filename.str());
  ShapeletObject *a = new ShapeletObject(average,averageBeta,xcentroid);
  a->setCartesianCoeffErrors(std_mean);
  
  filename.str("");
  filename << path << "average.sif";
  a->save(filename.str());
  delete a;
}

/// Creates shapelet image base on the dominant shapelet states
/// and their scatter, derived in Kelly & McKay (2003): astro-ph/0307395.
/// \f$N\f$ denotes the number of created galaxies.
void createShapeletImagesPCA(double betamin, double betamax,std::string path, int N) {
  NumMatrix<double> coeffs(9,9), sigma(9,9);
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
  createShapeletImages(coeffs,sigma,betamin,betamax,path,N);
}

void createLensedShapeletImages(double dx, NumMatrix<double>& kappa, NumMatrix<Complex>& gamma, NumMatrix<Complex>& F, NumMatrix<Complex>& G, std::string listfilename, std::string writeDirectory, int NOBJ) {
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
  NumMatrix<double> dGamma(2,2);
  for (int i = 0; i < N; i++) {
    for (int j =0 ; j < J; j++) {
      for (int n = 0; n < NOBJ; n++) {
	// random number from filelist
	randomnumber = (int)floor((double)filenumber*rand()/(RAND_MAX+1.0));
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
  if (argc != 5) {
    std::cout << "usage: createShapeletGalaxies <path> <beta_min> <beta_max> <NOBJ>" << std::endl;
    std::terminate();
  } 

  int NGLX = atoi(argv[4]);
  double betamin = atof(argv[2]), betamax = atof(argv[3]);
  createShapeletImagesPCA(betamin,betamax,argv[1],NGLX);
}
