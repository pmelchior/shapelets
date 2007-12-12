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
  TCLAP::CmdLine cmd("Decompose FITS file into shapelets", ' ', "0.2");
  TCLAP::UnlabeledValueArg<std::string> fileArg("file","FITS file to analyze",true,"","string", cmd);
  TCLAP::ValueArg<std::string> sifArg("p","prefix","Prefix of SIF files",true,"","string", cmd);
  TCLAP::ValueArg<std::string> confArg("c","config","ShapeLens configuration file",false,"","string", cmd);
  TCLAP::ValueArg<std::string> catArg("C","catalog","Name of catalog file (saved as ASCII file)",false,"","string",cmd);
  TCLAP::ValueArg<std::string> segmapArg("s","segmap","Name of segmentation map (stored as FITS file)",false,"","string", cmd);
  cmd.parse( argc, argv );
  
  // if a configuration file is provide, use it
  if (confArg.isSet())
    ShapeLensConfig sc(confArg.getValue());
  
  // open a FITS file and read given extension.
  // according to the cfitsio conventions, for th pHDU you can just
  // use 'filename', for extension you have to give 'filename[extension]'
  Frame* f = new Frame(fileArg.getValue());
  // measure noise background and subtract it by iterative sigma-clipping
  f->subtractBackground();

  // find objects within the image
  // the selection criteria are:
  // at least 50 pixels beyond 1.5 sigma above noise
  // and at least one pixel beyond 5 sigma (detection).
  f->findObjects();
  // return number of objects found
  unsigned int nobjects = f->getNumberOfObjects();

  // write segmentation map to FITS file if required
  // the number in this file are the running numbers for the objects
  // extracted below.
  if (segmapArg.isSet())
    f->getSegmentationMap().save(segmapArg.getValue());    
    
  // run through all objects
  // every file generated will have the appendix "_n", with n being the
  // object's running number
  for(int n = 1; n <= nobjects; n++) {

    // choose the actual object in the frame
    Object* obj = new Object(n);
    // "cut out" the object from whole frame and put it into Object obj
    f->fillObject(*obj);
    // dismiss objects with flag > 3 because of serious trouble
    // during the detection/segmentation process
    if (obj->getDetectionFlag() <= 3) {
      // create FITS file with the original pixel data of the object...
      std::ostringstream newname;
      newname << sifArg.getValue() << "_" << n << ".fits";
      std::string fitsname = newname.str();
      obj->save(fitsname);

      // actual decomposition is done here
      ShapeletObject sobj (*obj);
      // save results
      newname.str("");
      newname << sifArg.getValue() << "_" << n << ".sif";
      sobj.save(newname.str());

      // ... add shapelet model ...
      addFITSExtension(fitsname,"MODEL",obj->getGrid(),sobj.getModel());
      appendFITSHistory(fitsname+"[MODEL]",sobj.getHistory());
      // ... and residuals (data - model).
      NumVector<data_t> residuals = obj->getData();
      residuals -= sobj.getModel();
      addFITSExtension(fitsname,"RESIDUAL",obj->getGrid(),residuals);
    }
    delete obj;
  }
  if (catArg.isSet())
    f->getCatalog().save(catArg.getValue());

  delete f;
}
