#include <ShapeLens.h>
#include <tclap/CmdLine.h>

int main(int argc, char *argv[]) {
  TCLAP::CmdLine cmd("Create a FITS file from a shapelet model specified by SIF file", ' ', "0.2");
  TCLAP::ValueArg<data_t> rms("r","rms","RMS of pixel noise",false,0,"data_t", cmd);
  TCLAP::ValueArg<data_t> mean("m","mean","Mean of pixel noise",false,0,"data_t", cmd);
  TCLAP::SwitchArg poisson("p","poissonian","Poissonian noise model", cmd, false);
  TCLAP::SwitchArg history("H","history","Show hisory of SIF file", cmd, false);
  TCLAP::SwitchArg coeffs("c","coeffs","Save coefficients", cmd, false);
  TCLAP::SwitchArg errors("e","errors","Save coefficient errors", cmd, false);
  TCLAP::UnlabeledValueArg<std::string> input("input","Name of SIF file",true,"","string", cmd);
  TCLAP::UnlabeledValueArg<std::string> output("output","Name of FITS file",true,"","string", cmd);
  cmd.parse(argc,argv);

  ShapeletObject sobj(input.getValue());
  const Grid& grid = sobj.getGrid();
  NumVector<data_t> data = sobj.getModel();
  
  // add noise
  if (rms.isSet() || mean.isSet()) {
    if (poisson.isSet())
      addPoissonianNoise(data,mean.getValue());
    else
      addGaussianNoise(data,mean.getValue(),rms.getValue());
  }

  fitsfile* fptr = createFITSFile(output.getValue());
  writeFITSImage(fptr,grid,data,"MODEL");
  // append history
  if (history.isSet())
    appendFITSHistory(fptr,sobj.getHistory());
  // save coeffs and errors...
  if (coeffs.isSet())
    writeFITSImage(fptr,sobj.getCoeffs(),"COEFFS");
  if (errors.isSet())
    writeFITSImage(fptr,sobj.getDecompositionErrors(),"ERRORS");
  closeFITSFile(fptr);
}
  
							
