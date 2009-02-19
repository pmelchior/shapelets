#include <ShapeLens.h>
#include <tclap/CmdLine.h>
#include <fstream>

int main(int argc, char *argv[]) {
  TCLAP::CmdLine cmd("Perform shapelet decomposition of input file, save results and output shape catalog", ' ', "0.3");
  TCLAP::ValueArg<std::string> input("f","file","FITS file to analyze",true,"","string", cmd);
  TCLAP::ValueArg<std::string> prefix("p","prefix","Prefix of SIF files",true,"","string", cmd);
  TCLAP::ValueArg<std::string> config("c","config","ShapeLens configuration file",false,"","string", cmd);
  TCLAP::ValueArg<std::string> catalog("C","catalog","SExtractor catalog file",true,"","string", cmd);
  TCLAP::ValueArg<std::string> segmap("s","segmap","Name of segmentation map FITS file",true,"","string", cmd);
  TCLAP::ValueArg<std::string> weightmap("w","weightmap","Name of weight map FITS file",false,"","string", cmd);
  TCLAP::ValueArg<std::string> list("l","list","Name of file which lists the saved SIF files",false,"","string", cmd);
  TCLAP::SwitchArg model("m","model","Store object, model and residuals", cmd, false);
  cmd.parse( argc, argv );
  
  // set config options when given
  if (config.isSet())
    ShapeLensConfig sc(config.getValue());
  
  // open fits file to be analyzed
  SExFrame* f;
  if (weightmap.isSet())
    f = new SExFrame(input.getValue(),weightmap.getValue(),segmap.getValue(),catalog.getValue());
  else
    f = new SExFrame(input.getValue(),segmap.getValue(),catalog.getValue());
  f->subtractBackground();

  // if required: save a file which lists all stored SIF files
  std::ofstream listfile;
  if (list.isSet())
    listfile.open(list.getValue().c_str(),std::ios::out);
  
  // run through all objects in the catalog
  // this method ensures that only objects in the catalog are used
  // every file generated will have the appendix "_id", with id being the
  // object's running number
  const Catalog& cat = f->getCatalog();
  Catalog::const_iterator iter;
  for(iter = cat.begin(); iter != cat.end(); iter++) {
    // for clearity:
    unsigned long id = (*iter).first;
    // choose the actual object in the frame
    Object obj(id);
    // cut out object in place it in obj
    f->fillObject(obj);
    // dismiss objects with flags[i] = 1 for i >= 3 because of serious trouble
    // during the detection/segmentation process
    if ((*iter).second.FLAGS < 8) {

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
  // clean up
  if (list.isSet())
    listfile.close();
  delete f;
}
  
