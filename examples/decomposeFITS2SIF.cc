// the central header file for shapelet related work.
// you will need this almost always when working with the 'shapelets' library.
#include <ShapeletObject.h>
// the header file for the 'frame' library.
// since the 'frame' library is modular, this can be replaced by a file
// more appropriate to your needs.
#include <Frame.h>


// decomposeFITS2SIF:
// * open FITS file
// * estimate background noise mean and variance,
// * subtract noise
// * detect all objects large and bright enough to be significant
// * write segmentation map to FITS file
// * reject all objects truncated at the border (or otherwise troublesome)
// * decompose the remaining objects into shapelet states
// * store necessary information in SIF file
// * create FITS file with original object, shapelet model and residuals

int main(int argc, char *argv[]) {
  // read in necessary information from command line
  // complain if they are not provided
  if (!(argc == 3 || argc == 4) && !(argc == 7 || argc == 8)) {
    std::cout << "usage:" << std::endl;
    std::cout << "decomposeFITS2SIF <fits file[extension]> <sif file (without .sif)> [<regularization limit>]" << std::endl;
    std::cout << "decomposeFITS2SIF <fits file[extension]> <sif file (without .sif)> <nmaxLow> <nmaxHigh> <betaLow> <betaHigh> [<regularization limit>]" << std::endl;
  }
  else {
    std::string filename = argv[1];
    std::string sifprefix = argv[2];

    // if you want to see what the is actually going on, set this to 1.
    // if you want to have it quiet, set it to 0 (default).
    History::setVerbosity(1);

    // open a FITS file and read given extension.
    // according to the cfitsio conventions, for th pHDU you can just
    // use 'filename', for extension you have to give 'filename[extension]'
    Frame* f = new Frame(filename);
    // measure noise background and subtract it by iterative sigma-clipping
    f->subtractBackground();

    // find objects within the image
    // the selection criteria are:
    // at least 50 pixels beyond 1.5 sigma above noise
    // and at least one pixel beyond 5 sigma (detection).
    f->findObjects(50,1.5,5);
    // return number of objects found
    unsigned int nobjects = f->getNumberOfObjects();

    // write segmentation map to FITS file
    // the number in this file are the running numbers for the objects
    // extracted below.
    std::ostringstream newname;
    newname << sifprefix << "_seg.fits";
    // create set of keywords written to the FITS header (not mandatory).
    std::map<std::string, std::string> keys;
    keys["INFO"] = "segmentation map of " + filename;
    writeFITSFile(newname.str(),f->getGrid(),f->getObjectMap(),keys);
    
    
    // run through all objects
    // every file generated will have the appendix "_n", with n being the
    // object's running number
    for(int n = 1; n <= nobjects; n++) {
      std::ostringstream siffile;
      siffile << sifprefix << "_" << n << ".sif";

      // choose the actual object in the frame
      Object* obj = new Object(n);
      // "cut out" the object from whole frame and put it into Object obj
      f->fillObject(*obj);
      // employ gaussian statistics
      // this is faster and OK for faint objects which are background dominated
      // for bright ones, use "POISSONIAN".
      obj->setNoiseModel("GAUSSIAN");

      // dismiss objects with flag > 3 because of serious trouble
      // during the detection/segmentation process
      if (obj->getDetectionFlag() <= 3) {
	ShapeletObject *sobj;

	// automatic optimization: no bounds on beta and nmax
	if (argc == 3 || argc == 4) {
	  // when regularization should be done:
	  // specify the wanted limit on R
	  if (argc == 4){
	    ShapeletObject::DEFAULTS::REGULARIZE = 1;
	    ShapeletObject::DEFAULTS::REG_LIMIT = atof(argv[3]);
	  }
	  sobj =  new ShapeletObject(*obj);
	  sobj->save(siffile.str());
	}
	// use specific bounds for beta and nmax
	else if (argc == 7 || argc == 8) {
	  ShapeletObject::DEFAULTS::NMAX_LOW = atoi(argv[3]);
	  ShapeletObject::DEFAULTS::NMAX_HIGH = atoi(argv[4]);
	  ShapeletObject::DEFAULTS::BETA_LOW = atof(argv[5]);
	  ShapeletObject::DEFAULTS::BETA_HIGH = atof(argv[6]);
	  // when regularization should be done:
	  // specify the wanted limit on R
	  if (argc == 8){
	    ShapeletObject::DEFAULTS::REGULARIZE = 1;
	    ShapeletObject::DEFAULTS::REG_LIMIT = atof(argv[7]);
	  }
	  
	  // here comes the actual decomposition process.
	  // when this method is completed you can work entirely in
	  // shapelet space
	  sobj = new ShapeletObject(*obj);
	  // save all necessary information (shapelet parameters,
	  // coefficients, position of the object etc.) to a binary
	  // SIF file.
	  sobj->save(siffile.str());
	}

	// create FITS file with the original pixel data of the object...
	newname.str("");
	newname << sifprefix << "_" << n << ".fits";
	std::ostringstream info;
	info <<  "object " << n << " from " << filename;
	keys["INFO"] = info.str();
	writeFITSFile(newname.str(),obj->getGrid(),obj->getData(),keys);
	// ... add shapelet model ...
	addFITSExtension(newname.str(),"MODEL",obj->getGrid(),sobj->getModel());
 	// ... and residuals (data - model).
 	NumVector<double> residuals = obj->getData();
 	residuals -= sobj->getModel();
 	addFITSExtension(newname.str(),"RESIDUAL",obj->getGrid(),residuals);
	
	// clean up.
	delete sobj;
      }
      delete obj;
    }
    delete f;
  }
}
  
							
