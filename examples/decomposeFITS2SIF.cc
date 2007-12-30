// the central header file for shapelet related work.
// you will need this almost always when working with the 'shapelens' library.
#include <ShapeLens.h>
#include <tclap/CmdLine.h>

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
  TCLAP::CmdLine cmd("Decompose FITS file into shapelets", ' ', "0.3");
  TCLAP::UnlabeledValueArg<std::string> input("file","FITS file to analyze",true,"","string", cmd);
  TCLAP::ValueArg<std::string> prefix("p","prefix","Prefix of SIF files",true,"","string", cmd);
  TCLAP::ValueArg<std::string> config("c","config","ShapeLens configuration file",false,"","string", cmd);
  TCLAP::ValueArg<std::string> list("l","list","Name of file which lists the saved SIF files",false,"","string", cmd);
  TCLAP::ValueArg<std::string> catalog("C","catalog","Name of catalog file (saved as ASCII file)",false,"","string",cmd);
  TCLAP::ValueArg<std::string> segmap("s","segmap","Name of segmentation map (stored as FITS file)",false,"","string", cmd);
  TCLAP::SwitchArg model("m","model","Store object, model and residuals", cmd, false);
  cmd.parse( argc, argv );
  
  // if a configuration file is provide, use it
  if (config.isSet())
    ShapeLensConfig sc(config.getValue());
  
  // open a FITS file and read given extension.
  // according to the cfitsio conventions, for th pHDU you can just
  // use 'filename', for extension you have to give 'filename[extension]'
  Frame f(input.getValue());
  // measure noise background and subtract it by iterative sigma-clipping
  f.subtractBackground();

  // find objects within the image
  // the selection criteria are:
  // at least 50 pixels beyond 1.5 sigma above noise
  // and at least one pixel beyond 5 sigma (detection).
  f.findObjects();
  // return number of objects found
  unsigned int nobjects = f.getNumberOfObjects();

  // write segmentation map to FITS file if required
  // the number in this file are the running numbers for the objects
  // extracted below.
  if (segmap.isSet())
    f.getSegmentationMap().save(segmap.getValue());    
    
  // if required: save a file which lists all stored SIF files
  std::ofstream listfile;
  if (list.isSet())
    listfile.open(list.getValue().c_str(),std::ios::out);
  
  // run through all objects
  // every file generated will have the appendix "_n", with n being the
  // object's running number
  for(int n = 1; n <= nobjects; n++) {

    // choose the actual object in the frame
    Object obj(n);
    // "cut out" the object from whole frame and put it into Object obj
    f.fillObject(obj);
    // dismiss objects with flags[i] = 1 for i >= 3 because of serious trouble
    // during the detection/segmentation process
    if (obj.getDetectionFlags().to_ulong() < 8) {

      // actual decomposition is done here
      ShapeletObject sobj (obj);
      // save results
      std::ostringstream newname;
      newname << prefix.getValue() << "_" << n << ".sif";
      sobj.save(newname.str());
      if (list.isSet())
	listfile << newname.str() << std::endl;

      // if a model should be stored, do it here:
      if (model.isSet()) {
	// create FITS file with the original pixel data of the object...
	newname.str("");
	newname << prefix.getValue() << "_" << n << ".fits";
	std::string fitsname = newname.str();
	obj.save(fitsname);
	// ... add shapelet model ...
	fitsfile* fptr = openFITSFile(fitsname,1);
	writeFITSImage(fptr,obj.getGrid(),sobj.getModel(),"MODEL");
	appendFITSHistory(fptr,sobj.getHistory());
	// ... and residuals (data - model).
	NumVector<data_t> residuals = obj.getData();
	residuals -= sobj.getModel();
	writeFITSImage(fptr,obj.getGrid(),residuals,"RESIDUAL");
	closeFITSFile(fptr);
      }
    }
  }
  if (catalog.isSet())
    f.getCatalog().save(catalog.getValue());

  if (list.isSet())
    listfile.close();
}
