#include <ShapeLens.h>
#include <tclap/CmdLine.h>

// convertSIF2PPM:
// create a PPM file from a shapelet model specified by SIF file.

int main(int argc, char *argv[]) {
  TCLAP::CmdLine cmd("Create a PPM file from a shapelet model specified by SIF file", ' ', "0.1");
  TCLAP::ValueArg<data_t> rms("r","rms","RMS of pixel noise",false,0,"data_t", cmd);
  TCLAP::ValueArg<data_t> mean("m","mean","Mean of pixel noise",false,0,"data_t", cmd);
  TCLAP::SwitchArg poisson("p","poissonian","Poissonian noise model", cmd, false);
  TCLAP::UnlabeledValueArg<std::string> input("input","Name of SIF file",true,"","string", cmd);
  TCLAP::UnlabeledValueArg<std::string> output("output","Name of PPM file",true,"","string", cmd);
  cmd.parse(argc,argv);

  // open sif file
  ShapeletObject sobj(input.getValue());
  // construct shapelet model from coeffs stored in sobj
  Image<data_t> data = sobj.getModel();
  // add noise if required

  // add noise
  if (rms.isSet() || mean.isSet()) {
    if (poisson.isSet())
      IO::addPoissonianNoise(data,mean.getValue());
    else
      IO::addGaussianNoise(data,mean.getValue(),rms.getValue());
  }

  // write model to ppm file:
  // colormodel in {"RED", "BLUE", "GREEN", "GRAY", "WARM"}
  // data scaling in {"LINEAR", "SQUARE_ROOT", "LOGARITHMIC"} 
  // for the data range between min() and max()
  IO::writePPMImage(output.getValue(),"SPECTRUM","LINEAR",data.min(),data.max(),data.grid,data);
}
  
							
