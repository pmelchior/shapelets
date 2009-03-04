#include <ShapeLens.h>
#include <tclap/CmdLine.h>

using namespace shapelens;

// printSIFFile:
// open SIF file and print header (+ decomposition history + 
// shapelet coefficients/errors) to stdout.

int main (int argc, char *argv[]) {
  TCLAP::CmdLine cmd("Print SIF file information", ' ', "0.2");
  TCLAP::SwitchArg historySwitch("H","history","Show hisory of SIF file", cmd, false);
  TCLAP::SwitchArg coefficientSwitch("c","coefficients","Show coefficients of SIF file", cmd, false);
  TCLAP::SwitchArg polarSwitch("p","polar_coeffs","Show polar coefficients of SIF file", cmd, false);
  TCLAP::SwitchArg errorSwitch("e","errors","Show coefficient errors of SIF file", cmd, false);
  TCLAP::SwitchArg covSwitch("C","covariance","Show covariance matrix of SIF file", cmd, false);
  TCLAP::SwitchArg sigSwitch("s","significance","Show significance (= S/N) of the coefficients of SIF file", cmd, false);
  TCLAP::UnlabeledValueArg<std::string> nameArg("name","Name of SIF file",true,"","string", cmd);
  cmd.parse(argc,argv);

  // open ShapeletObject from SIFFile
  ShapeletObject sobj(nameArg.getValue());
  // always print out general infos:
  std::cout << "SIF header information" << std::endl;
  std::cout << "Filename:\t" << sobj.getBaseFilename() << std::endl;
  std::cout << "Coefficients:\t" << sobj.getCoeffs().size() << " (n_max = " << sobj.getNMax() << ")" << std::endl;
  std::cout << "Coff. errors:\t" << sobj.getCovarianceMatrix().getRows() << " x " << sobj.getCovarianceMatrix().getColumns() << std::endl;;
  std::cout << "Beta:\t\t" << sobj.getBeta() << std::endl;
  std::cout << "Centroid:\t" << (sobj.getCentroid())(0) << "/" << (sobj.getCentroid())(1) << std::endl;
  std::cout << "Grid:\t\t" << (sobj.getGrid()).getStartPosition(0) << ".." << (sobj.getGrid()).getStopPosition(0) << ", " << (sobj.getGrid()).getStartPosition(1) << ".." << (sobj.getGrid()).getStopPosition(1) << std::endl;
  std::cout << "Chi^2:\t\t" << sobj.getChiSquare() << std::endl;
  std::cout << "Flags:\t\t" << sobj.getFlags().to_string<char,std::char_traits<char>,std::allocator<char> >().insert(8," ") << std::endl;
  std::cout << std::endl;

  // if required, print out decomposition history
  if (historySwitch.isSet())
    std::cout << "History:" << std::endl << sobj.getHistory() << std::endl;
  // if required, print out shapelet coeffs (+ errors)
  if (coefficientSwitch.isSet())
    std::cout << "Cartesian coefficients:" << std::endl << sobj.getCoeffs() << std::endl << std::endl;
  if (polarSwitch.isSet())
    std::cout << "Polar coefficients:" << std::endl << sobj.getPolarCoeffs() << std::endl << std::endl;
  if (covSwitch.isSet())
    std::cout << "Cartesian coefficient covariance:" << std::endl << sobj.getCovarianceMatrix() << std::endl << std::endl;
  if (errorSwitch.isSet())
    std::cout <<  "Cartesian coefficient errors:" << std::endl << sobj.getErrors() << std::endl << std::endl;
  if (sigSwitch.isSet()) {
    CoefficientVector<data_t> sig = sobj.getErrors();
    const CoefficientVector<data_t>& coeffs = sobj.getCoeffs();
    for (unsigned int i=0; i< sig.size(); i++)
      sig(i) = fabs(coeffs(i))/sig(i);
    std::cout <<  "Cartesian coefficient significance:" << std::endl << sig << std::endl;
    std::cout <<  "mean significance: " << sig.mean() << std::endl << std::endl;
  }
}

