#include <ShapeLens.h>
#include <tclap/CmdLine.h>

// convertSIF2PPM:
// create a PPM file from a shapelet model specified by SIF file.

int main(int argc, char *argv[]) {
  TCLAP::CmdLine cmd("Create a PPM file from a shapelet model specified by SIF file", ' ', "0.1");
  TCLAP::ValueArg<data_t> rmsArg("r","rms","RMS of pixel noise",false,0,"data_t", cmd);
  TCLAP::ValueArg<data_t> meanArg("m","mean","Mean of pixel noise",false,0,"data_t", cmd);
  TCLAP::SwitchArg poissonSwitch("p","poissonian","Poissonian noise model", cmd, false);
  TCLAP::UnlabeledValueArg<std::string> inArg("input","Name of SIF file",true,"","string", cmd);
  TCLAP::UnlabeledValueArg<std::string> outArg("output","Name of PPM file",true,"","string", cmd);
  cmd.parse(argc,argv);

  // open sif file
  ShapeletObject sobj(inArg.getValue());
  // get grid from sobj
  const Grid& grid = sobj.getGrid();
  // construct shapelet model from coeffs stored in sobj
  NumVector<data_t> data = sobj.getModel();
  // add noise if required
  if (poissonSwitch.isSet())
    addPoissonianNoise(data,meanArg.getValue(),rmsArg.getValue());
  else
    addGaussianNoise(data,meanArg.getValue(),rmsArg.getValue());

  // write model to ppm file:
  // colormodel in {"RED", "BLUE", "GREEN", "GRAY", "WARM"}
  // data scaling in {"LINEAR", "SQUARE_ROOT", "LOGARITHMIC"} 
  // for the data range between min() and max()
  writePPMImage(outArg.getValue(),"SPECTRUM","LINEAR",data.min(),data.max(),grid,data);
}
  
							
