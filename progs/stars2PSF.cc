#include <shapelens/ShapeLens.h>
#include <shapelens/lensing/DEIMOS.h>
#include <tclap/CmdLine.h>

using namespace shapelens;

data_t xi(data_t r) {
  return exp(-r/2100);
}

int main(int argc, char* argv[]) {
  TCLAP::CmdLine cmd("Measure PSF moments from stars and interpolate in between", ' ', "0.2");
  TCLAP::ValueArg<std::string> cat("c","catalog","Catalog file", true, "","string",cmd);
  TCLAP::ValueArg<std::string> file("f","file","Image file", true, "","string",cmd);
  TCLAP::ValueArg<int> C("C","correction_order","DEIMOS correction order", true, 4,"int",cmd);
  TCLAP::ValueArg<data_t> s("s","scale","DEIMOS weighting scale", true, 3.,"data_t",cmd);
  TCLAP::SwitchArg flexed("F","flexed","Enable flexion in DEIMOS", cmd, false);
  TCLAP::SwitchArg printMoments("m","print_moments","Print moments", cmd, false);
  TCLAP::SwitchArg printErrors("e","print_errors","Print moment errors", cmd, false);
  TCLAP::ValueArg<std::string> weight("w","weightmap","Weight map file", false, "","string");
  TCLAP::ValueArg<data_t> noise("n","noise_rms","RMS of the background noise", true, 0.,"data_t");
  cmd.xorAdd(weight, noise);
  TCLAP::ValueArg<std::string> saveObj("S","save_objects","Save objects in this file", false, "", "string",cmd);
  TCLAP::ValueArg<std::string> average("a","average","Save average moments to file", false, "", "string",cmd);
  cmd.parse(argc,argv);

  // shear/flexion?
  int N = 2;
  if (flexed.getValue())
    N = 4;


  // open file pointers and catalog
  HugeFrame* frame;
  if (weight.isSet())
    frame = new HugeFrame(file.getValue(), weight.getValue(), cat.getValue());
  else
    frame = new HugeFrame(file.getValue(), cat.getValue());
  frame->setNoiseMeanRMS(0,noise.getValue());

  Object obj;
  fitsfile* fptr;
  if (saveObj.isSet())
    fptr = IO::createFITSFile(saveObj.getValue());

  DEIMOS avg;

  for (Catalog::const_iterator iter = frame->getCatalog().begin(); 
       iter != frame->getCatalog().end();
       iter++) {
    frame->fillObject(obj,iter);
    DEIMOS d(obj, N, C.getValue(), s.getValue(), flexed.getValue());
    
    if (iter == frame->getCatalog().begin())
      avg = d;
    else {
      avg.mo += d.mo; // flux weighted average
      avg.S += d.S;
      avg.eps += d.eps;
      if (flexed.getValue())
	avg.G += d.G;
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
    if (flexed.getValue()) {
      data_t trQ = d.mo(2,0) + d.mo(0,2);
      data_t xi = d.mo(4,0) + 2*d.mo(2,2) + d.mo(0,4);
      data_t mu = trQ*trQ/xi;
      std::cout << "\t" << real(d.delta()) << "\t" << imag(d.delta());
      std::cout << "\t" << mu;
    }
    std::cout << std::endl;
    if (saveObj.isSet())
      IO::writeFITSImage(fptr,obj);
  }

  // clean up
  delete frame;
  if (saveObj.isSet())
    IO::closeFITSFile(fptr);

  if (average.isSet()) {
    // flux-weighted average
    avg.mo /= frame->getCatalog().size();
    avg.S /= frame->getCatalog().size();
    // flux normalization
    data_t flux = avg.mo(0,0);
    avg.mo /= flux;
    avg.S /= flux*flux;
    // plain average
    avg.eps /= frame->getCatalog().size();
    avg.G /= frame->getCatalog().size();
    avg.id = 0;
    avg.save(average.getValue());
  }
    

}
