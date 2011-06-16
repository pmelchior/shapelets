#include <shapelens/ShapeLens.h>
#include <shapelens/lensing/DEIMOS.h>
#include <tclap/CmdLine.h>

using namespace shapelens;

int main(int argc, char* argv[]) {
  TCLAP::CmdLine cmd("Measure galactic moments", ' ', "0.3");
  TCLAP::ValueArg<std::string> cat("c","catalog","Catalog file", true, "","string",cmd);
  TCLAP::ValueArg<std::string> weight("w","weightmap","Weight map file", false, "","string");
  TCLAP::ValueArg<data_t> noise("n","noise_rms","RMS of the background noise", true, 0.,"data_t");
  cmd.xorAdd(weight, noise);
  TCLAP::ValueArg<std::string> psffile("p","psf_file","PSF moments file", true, "","string",cmd);
  TCLAP::ValueArg<std::string> file("f","file","Image file", true, "","string",cmd);
  TCLAP::ValueArg<int> order("N","moment_order","Moment order", false, 2,"int",cmd);
  TCLAP::ValueArg<int> C("C","correction_order","DEIMOS correction order", true, 4,"int",cmd);
  TCLAP::ValueArg<data_t> s("s","scale","DEIMOS weighting scale", true, 3.,"data_t",cmd);
  TCLAP::SwitchArg flexed("F","flexed","Enable flexion in DEIMOS", cmd, false);
  TCLAP::SwitchArg printMoments("m","print_moments","Print moments", cmd, false);
  TCLAP::SwitchArg printErrors("e","print_errors","Print moment errors", cmd, false);
  TCLAP::SwitchArg usewcs("u","use_wcs","Use WCS from FITS header", cmd, false);
  cmd.parse(argc,argv);

  // shear/flexion?
  int N = order.getValue();
  if (flexed.getValue())
    N = std::min(4,N);

  // set WCS if requested
  ShapeLensConfig::USE_WCS = usewcs.getValue();

  // open file pointers and catalog
  HugeFrame* frame;
  if (weight.isSet())
    frame = new HugeFrame(file.getValue(), weight.getValue(), cat.getValue());
  else
    frame = new HugeFrame(file.getValue(), cat.getValue());
  frame->setNoiseMeanRMS(0,noise.getValue());

  // open PSF moment maps
  // only moments higher than 1 are stored
  // assumption: mom_00 = 1, mom_01 = mom_10 = 0

  // DEIMOS psf;
//   psf.mo.setOrder(N);
//   psf.mo(0,0) = 1;
//   psf.mo(0,1) = psf.mo(1,0) = 0;
//   fitsfile* fptr_p = IO::openFITSFile(psffile.getValue());
//   int B;
//   IO::readFITSKeyword(fptr_p,"B",B); // kriging box size
//   int m = 3;
//   if (N == 4)
//     m = 12;
//   Image<data_t> map;
//   std::vector<Image<data_t> > mom_psf;
//   for (int i = 0; i < m; i++) {
//     IO::moveToFITSExtension(fptr_p,i+1);
//     IO::readFITSImage(fptr_p,map);
//     mom_psf.push_back(map);
//   }
//   IO::closeFITSFile(fptr_p);
  DEIMOS p(psffile.getValue());
  DEIMOS::PSFMultiScale psf;
  psf.insert(s.getValue(), p.mo);

  Object obj;
  for (Catalog::const_iterator iter = frame->getCatalog().begin(); 
       iter != frame->getCatalog().end();
       iter++) {
    frame->fillObject(obj,iter);
    DEIMOS d(obj, psf, N, C.getValue(), s.getValue(), flexed.getValue());

    // // fill psf moments from moment maps
//     int pos1 = int(floor(obj.centroid(0)))/B;
//     int pos2 = int(floor(obj.centroid(1)))/B;
//     for (int i=3; i < psf.mo.size(); i++)
//       psf.mo(i) = (mom_psf[i-3])(pos1,pos2);

    std::cout << iter->first << "\t" << obj.centroid(0) << "\t" << obj.centroid(1) << "\t" ;
    if (printMoments.getValue()) 
      for (int i=0; i < d.mo.size(); i++)
	std::cout << d.mo(i) << "\t";
    if (printErrors.getValue()) {
      Moments mo_noise = d.getMomentErrors();
      for (int i=0; i < d.mo.size(); i++)
	std::cout << sqrt(mo_noise(i)) << "\t";
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
  }

  // clean up
  delete frame;

}
