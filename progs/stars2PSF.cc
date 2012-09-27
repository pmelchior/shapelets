#include <shapelens/ShapeLens.h>
#include <shapelens/lensing/DEIMOSElliptical.h>
#include <tclap/CmdLine.h>

using namespace shapelens;

int main(int argc, char* argv[]) {
  TCLAP::CmdLine cmd("Measure PSF moments from stars and interpolate in between", ' ', "0.4");
  TCLAP::ValueArg<std::string> cat("c","catalog","Catalog file", true, "","string",cmd);
  TCLAP::ValueArg<std::string> segmap("S","segmap","Segmentation maps", false, "", "string",cmd);
  TCLAP::ValueArg<std::string> weight("w","weightmap","Weight map file", false, "","string");
  TCLAP::ValueArg<std::string> file("f","file","Image file", true, "","string",cmd);
  TCLAP::ValueArg<int> order("N","moment_order","Moment order", false, 2,"int",cmd);
  TCLAP::ValueArg<int> C("C","correction_order","DEIMOS correction order", true, 4,"int",cmd);
  TCLAP::ValueArg<data_t> s("s","scale","DEIMOS weighting scale", true, 3.,"data_t",cmd);
  TCLAP::SwitchArg printMoments("m","print_moments","Print moments", cmd, false);
  TCLAP::SwitchArg printErrors("e","print_errors","Print moment errors", cmd, false);
  TCLAP::ValueArg<data_t> noise("n","noise_rms","RMS of the background noise", false, 0.,"data_t");
  cmd.xorAdd(weight, noise);
  TCLAP::ValueArg<std::string> saveObj("o","save_objects","Save objects in this file", false, "", "string",cmd);
  TCLAP::ValueArg<std::string> average("a","average","Save average moments to file", false, "", "string",cmd);
  TCLAP::SwitchArg usewcs("u","use_wcs","Use WCS from FITS header", cmd, false);
  TCLAP::SwitchArg verbose("v","verbose","Verbose output", cmd, false);
  cmd.parse(argc,argv);

  int N = order.getValue();

  // set WCS if requested
  ShapeLensConfig::USE_WCS = usewcs.getValue();
  ShapeLensConfig::CHECK_OBJECT = true;
  ShapeLensConfig::VERBOSITY = verbose.getValue();

  // open file pointers and catalog
  SExFrame frame(file.getValue(), cat.getValue(), segmap.getValue(), weight.getValue());
  frame.setNoiseMeanRMS(0,noise.getValue());

  Object obj;
  fitsfile* fptr;
  if (saveObj.isSet())
    fptr = IO::createFITSFile(saveObj.getValue());

  DEIMOSElliptical avg(N);

  for (Catalog::const_iterator iter = frame.getCatalog().begin(); 
       iter != frame.getCatalog().end();
       iter++) {
    frame.fillObject(obj,iter);
    if (iter->first == 1)
      obj.save("first.fits");
    DEIMOSElliptical d(obj, N, C.getValue(), s.getValue());
    
    if (iter == frame.getCatalog().begin())
      avg = d;
    else {
      avg.mo += d.mo; // flux weighted average
      avg.S += d.S;
      avg.eps += d.eps;
      avg.scale += d.scale;
    }
    
    std::cout << iter->first << "\t" << obj.centroid(0) << "\t" << obj.centroid(1) << "\t" ;
    if (printMoments.getValue()) 
      for (int i=0; i < d.mo.size(); i++)
	std::cout << d.mo(i) << "\t";
    if (printErrors.getValue()) {
      Moments mo_noise = d.getMomentErrors();
      for (int i=0; i < d.mo.size(); i++)
	std::cout << mo_noise(i) << "\t";
    }
    std::cout << real(d.epsilon()) << "\t" << imag(d.epsilon());
    std::cout << "\t" << d.SN[s.getValue()] << "\t" << d.flags.to_string() << std::endl;
    if (saveObj.isSet())
      IO::writeFITSImage(fptr,obj);
  }
  
  // clean up
  if (saveObj.isSet())
    IO::closeFITSFile(fptr);

  if (average.isSet()) {
    // flux-weighted average
    avg.mo /= frame.getCatalog().size();
    avg.S /= frame.getCatalog().size();
    // flux normalization
    data_t flux = avg.mo(0,0);
    avg.mo /= flux;
    avg.S /= flux*flux;
    // plain average
    avg.eps /= frame.getCatalog().size();
    avg.scale /= frame.getCatalog().size();
    avg.save(average.getValue());
  }
    

}
