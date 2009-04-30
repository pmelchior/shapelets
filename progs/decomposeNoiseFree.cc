#include <ShapeLens.h>
#include <tclap/CmdLine.h>

using namespace shapelens;

int main(int argc, char *argv[]) {
  TCLAP::CmdLine cmd("Decompose FITS file with nosie-free data", ' ', "0.1");
  TCLAP::UnlabeledValueArg<std::string> input("file","FITS file to analyze",true,"","string", cmd);
  TCLAP::UnlabeledValueArg<std::string> output("output","Name of SIF file",true,"","string", cmd);
  TCLAP::SwitchArg model("m","model","Create FITS file with model and residuals", cmd, false);
  TCLAP::ValueArg<int> sidelength("L","sidelength","Sidelength of object's box",false,100,"int", cmd);
  TCLAP::ValueArg<data_t> noise_rms("r","noise_rms","RMS of noise for normalizing chi2",false,1e-10,"data_t",cmd);
  TCLAP::ValueArg<int> nmax("n","n_max","Maximum shapelet order",false,24,"int",cmd);
  TCLAP::ValueArg<std::string> config("c","config","ShapeLens configuration file",false,"","string", cmd);
  cmd.parse(argc,argv);
  
  // if a configuration file is provide, use it
  if (config.isSet())
    ShapeLensConfig sc(config.getValue());
  if (nmax.isSet())
    ShapeLensConfig::NMAX_LOW = ShapeLensConfig::NMAX_HIGH = nmax.getValue();

  // form object from image data
  Object obj;
  obj = Image<data_t>(input.getValue());
  obj.computeFlux();
  obj.computeCentroid();
  if (sidelength.isSet()) {
    int L = sidelength.getValue();
    Point2D<int> p1,p2;
    p1(0) = obj.centroid(0) - L/2;
    p1(1) = obj.centroid(1) - L/2;
    p2(0) = obj.centroid(0) + L/2;
    p2(1) = obj.centroid(1) + L/2;
    Image<data_t> sub;
    obj.slice(sub,p1,p2);
    obj = sub;
  }
  obj.noise_mean = 0;
  obj.noise_rms = noise_rms.getValue(); // just a number to obtain sensible chi2
  
  // decompose it
  ShapeletObject sobj(obj);
  
  // save sif file
  sobj.save(output.getValue());

  if (model.isSet()) {
    std::string filename = output.getValue();
    size_t l = filename.find_last_of(".");
    if (l != std::string::npos) {
      filename = filename.replace(l,5,".fits");
      fitsfile* fptr = IO::createFITSFile(filename);
      IO::writeFITSImage(fptr,obj);
      // ... add shapelet model ...
      IO::writeFITSImage(fptr,sobj.getModel(),"MODEL");
      IO::appendFITSHistory(fptr,sobj.getHistory());
      // ... and residuals (obj - model).
      obj -= sobj.getModel();
      IO::writeFITSImage(fptr,obj,"RESIDUAL");
      IO::closeFITSFile(fptr);
    }
  }
}
