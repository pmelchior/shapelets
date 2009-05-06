#ifdef HAS_FFTW3

#include <shapelens/utils/FFT.h>
#include <shapelens/utils/IO.h>
#include <shapelens/frame/Object.h>
#include <tclap/CmdLine.h> 

using namespace shapelens;

int main(int argc, char *argv[]) {
  TCLAP::CmdLine cmd("Convolve FitsObject with anoher FitsObject in Real Space", ' ', "0.1");
  TCLAP::ValueArg<std::string> input("i","input","input FITS file",true,"","std::string", cmd);
  TCLAP::ValueArg<std::string> kernel ("k","kernel","kernel FITS file",true,"","std::string",cmd);
  TCLAP::ValueArg<std::string> output ("o","output","output FITS file",true,"","std::string",cmd);
  TCLAP::ValueArg<data_t> rms("r","rms","RMS of pixel noise",false,0,"data_t", cmd);
  TCLAP::ValueArg<data_t> mean("m","mean","Mean of pixel noise",false,0,"data_t", cmd);
  TCLAP::SwitchArg poisson("p","poissonian","Poissonian noise model", cmd, false);
  cmd.parse(argc,argv);

  // Load the input FITS files for data and kernel:
  //Image<data_t> im_input, im_kernel;
  Object im_input, im_kernel;
  fitsfile* fptr = IO::openFITSFile(input.getValue());
  IO::readFITSImage(fptr,im_input);
  IO::closeFITSFile(fptr);
  fptr = IO::openFITSFile(kernel.getValue());
  IO::readFITSImage(fptr,im_kernel);
  IO::closeFITSFile(fptr);

  Object im_output = im_input;
  im_output.convolve(im_kernel);
	
  // add noise
  if (rms.isSet() || mean.isSet()) {
    if (poisson.isSet())
      IO::addPoissonianNoise(im_input,mean.getValue());
    else
      IO::addGaussianNoise(im_input,mean.getValue(),rms.getValue());
  }
	
  // Save the output matrix to the FITS file given by user:
  fptr = IO::createFITSFile(output.getValue());
  IO::writeFITSImage(fptr, im_input);
  IO::writeFITSImage(fptr, im_kernel);
  IO::writeFITSImage(fptr, im_output);
  IO::closeFITSFile(fptr);
}
#else
int main(){
  return 0;
}
#endif
