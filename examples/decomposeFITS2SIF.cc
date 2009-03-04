// the central header file for shapelet related work.
// you will need this almost always when working with the 'shapelens' library.
#include <ShapeLens.h>
#include <tclap/CmdLine.h>

using namespace shapelens;

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
  TCLAP::CmdLine cmd("Decompose FITS file into shapelets", ' ', "0.5");
  TCLAP::UnlabeledValueArg<std::string> input("file","FITS file to analyze",true,"","string", cmd);
  TCLAP::ValueArg<std::string> prefix("p","prefix","Prefix of SIF files",true,"","string", cmd);
  TCLAP::ValueArg<std::string> config("c","config","ShapeLens configuration file",false,"","string", cmd);
  TCLAP::ValueArg<std::string> list("l","list","Name of file which lists the saved SIF files",false,"","string", cmd);
  TCLAP::ValueArg<std::string> catalog("C","catalog","Name of catalog file (saved as ASCII file)",false,"","string",cmd);
  TCLAP::ValueArg<std::string> segmap("s","segmap","Name of segmentation map (stored as FITS file)",false,"","string", cmd);
  TCLAP::ValueArg<std::string> weightmap("w","weightmap","Name of weight map FITS file",false,"","string", cmd);
  TCLAP::SwitchArg model("m","model","Store object, model and residuals", cmd, false);
  cmd.parse( argc, argv );
  
  // if a configuration file is provide, use it
  if (config.isSet())
    ShapeLensConfig sc(config.getValue());
  
  // open a FITS file and read given extension.
  // according to the cfitsio conventions, for th pHDU you can just
  // use 'filename', for extension you have to give 'filename[extension]'
  Frame* f;
  if (weightmap.isSet())
    f = new Frame(input.getValue(),weightmap.getValue());
  else
    f = new Frame(input.getValue());

  // measure noise background and subtract it by iterative sigma-clipping
  f->subtractBackground();

  // find objects within the image
  // the selection criteria are:
  // at least 50 pixels beyond 1.5 sigma above noise
  // and at least one pixel beyond 5 sigma (detection).
  f->findObjects();

  // write segmentation map to FITS file if required
  // the number in this file are the running numbers for the objects
  // extracted below.
  if (segmap.isSet())
    f->getSegmentationMap().save(segmap.getValue());    
  
  // if noisemodel is COVARIANCE: measure pixel correlations
  CorrelationFunction xi;
  if (ShapeLensConfig::NOISEMODEL == "COVARIANCE")
    xi = f->computeCorrelationFunction(2);

  // if required: save a file which lists all stored SIF files
  std::ofstream listfile;
  if (list.isSet())
    listfile.open(list.getValue().c_str(),std::ios::out);
  
  // run through all objects
  // every file generated will have the appendix "_n", with n being the
  // object's running number
  const Catalog& cat = f->getCatalog();
  Catalog::const_iterator iter;
  for(iter = cat.begin(); iter != cat.end(); iter++) {
    // for clearity:
    unsigned long id = (*iter).first;
    // choose the actual object in the frame:
    // "cut out" the object from whole frame and put it into Object obj
    Object obj;
    f->fillObject(obj,iter);
    // set objects correlation to the one measure from f
    if (ShapeLensConfig::NOISEMODEL == "COVARIANCE")
      obj.xi = xi;
    // dismiss objects with a flag set from detection/segmentation process
    if ((*iter).second.FLAGS == 0) {

      // actual decomposition is done here
      ShapeletObject sobj (obj);
      // save results
      std::ostringstream newname;
      newname << prefix.getValue() << "_" << id << ".sif";
      sobj.save(newname.str());
      if (list.isSet())
	listfile << newname.str() << std::endl;

      // if a model should be stored, do it here:
      if (model.isSet()) {
	// create FITS file with the original pixel data of the object...
	newname.str("");
	newname << prefix.getValue() << "_" << id << ".fits";
	std::string fitsname = newname.str();
	obj.save(fitsname);
	// ... add shapelet model ...
	fitsfile* fptr = IO::openFITSFile(fitsname,1);
	IO::writeFITSImage(fptr,sobj.getModel(),"MODEL");
	IO::appendFITSHistory(fptr,sobj.getHistory());
	// ... and residuals (obj - model).
	obj -= sobj.getModel();
	IO::writeFITSImage(fptr,obj,"RESIDUAL");
	IO::closeFITSFile(fptr);
      }
    }
  }
  // save catalog when demanded
  if (catalog.isSet())
    cat.save(catalog.getValue());

  // clean up
  if (list.isSet())
    listfile.close();
  delete f;
}
