#include <ShapeLens.h>
#include <tclap/CmdLine.h>

// printSIFFile:
// open SIF file and print header (+ decomposition history + 
// shapelet coefficients/errors) to stdout.

int main (int argc, char *argv[]) {
  TCLAP::CmdLine cmd("Print SIF file information", ' ', "0.1");
  TCLAP::SwitchArg historySwitch("H","history","Show hisory of SIF file", cmd, false);
  TCLAP::SwitchArg coefficientSwitch("c","coefficients","Show coefficients of SIF file", cmd, false);
  TCLAP::SwitchArg errorSwitch("e","errors","Show coefficient errors of SIF file", cmd, false);
  TCLAP::UnlabeledValueArg<std::string> nameArg("name","Name of SIF file",true,"","string", cmd);
  cmd.parse(argc,argv);

  // open ShapeletObject from SIFFile
  ShapeletObject sobj(nameArg.getValue());
  // always print out general infos:
  std::cout << "SIF header information" << std::endl;
  std::cout << "Filename:\t" << sobj.getBaseFilename() << std::endl;
  std::cout << "Coefficients:\t" << sobj.getCartesianCoeffs().getRows() <<"x" <<  sobj.getCartesianCoeffs().getColumns() << std::endl;
  std::cout << "Coff. errors:\t" << sobj.getDecompositionErrors().getRows() <<"x" << sobj.getDecompositionErrors().getRows() << std::endl;;
  std::cout << "Beta:\t\t" << sobj.getBeta() << std::endl;
  std::cout << "Centroid:\t" << (sobj.getCentroid())(0) << "/" << (sobj.getCentroid())(1) << std::endl;
  std::cout << "Grid:\t\t" << (sobj.getGrid()).getStartPosition(0) << ".." << (sobj.getGrid()).getStopPosition(0) << ", " << (sobj.getGrid()).getStartPosition(1) << ".." << (sobj.getGrid()).getStopPosition(1) << std::endl;
  std::cout << "Chi^2:\t\t" << sobj.getDecompositionChiSquare() << std::endl;
  std::cout << "Flags:\t\t" << (sobj.getFlags()).to_ulong() << std::endl;
  std::cout << "Regularized:\t" << ShapeLensConfig::REGULARIZE;
  if (ShapeLensConfig::REGULARIZE)
    std::cout << " (R = " << sobj.getRegularizationR() << ")" << std::endl;
  else 
    std::cout << std::endl;

  // if required, print out decomposition history
  if (historySwitch.isSet())
    std::cout << sobj.getHistory();
  // if required, print out shapelet coeffs (+ errors)
  if (coefficientSwitch.isSet())
    std::cout << "Cartesian coefficients:" << std::endl << sobj.getCartesianCoeffs() << std::endl;
  if (errorSwitch.isSet())
    std::cout <<  "Cartesian coefficient errors:" << std::endl << sobj.getDecompositionErrors() << std::endl;
  
}

