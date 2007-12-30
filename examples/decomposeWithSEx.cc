#include <ShapeLens.h>
#include <tclap/CmdLine.h>
#include <fstream>

int main(int argc, char *argv[]) {
  TCLAP::CmdLine cmd("Perform shapelet decomposition of input file, save results and output shape catalog", ' ', "0.2");
  TCLAP::ValueArg<std::string> input("f","file","FITS file to analyze",true,"","string", cmd);
  TCLAP::ValueArg<std::string> prefix("p","prefix","Prefix of SIF files",true,"","string", cmd);
  TCLAP::ValueArg<std::string> config("c","config","ShapeLens configuration file",false,"","string", cmd);
  TCLAP::ValueArg<std::string> catalog("C","catalog","SExtractor catalog file",true,"","string", cmd);
  TCLAP::ValueArg<std::string> segmap("s","segmap","Name of segmentation map FITS file",true,"","string", cmd);
  TCLAP::ValueArg<std::string> list("l","list","Name of file which lists the saved SIF files",false,"","string", cmd);
  TCLAP::SwitchArg model("m","model","Store object, model and residuals", cmd, false);
  cmd.parse( argc, argv );
  
  // set config options when given
  if (config.isSet())
    ShapeLensConfig sc(config.getValue());
  
  // open fits file to be analyzed
  SExFrame f(input.getValue());
  f.readCatalog(catalog.getValue());
  f.readSegmentationMap(segmap.getValue());
  f.subtractBackground();

  // if required: save a file which lists all stored SIF files
  std::ofstream listfile;
  if (list.isSet())
    listfile.open(list.getValue().c_str(),std::ios::out);
  

  // iterate through detected objects
  for (int n = 1; n <= f.getNumberOfObjects(); n++) {
    // create empty object container with ID = n
    Object obj(n);
    // cut out object in place it in obj
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
  if (list.isSet())
    listfile.close();
}
  
