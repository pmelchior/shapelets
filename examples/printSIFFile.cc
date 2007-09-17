#include <ShapeLens.h>

// printSIFFile:
// open SIF file and print header (+ decomposition history + 
// shapelet coefficients/errors) to stdout.

int main (int argc, char *argv[]) {
  if (argc != 4)
    std::cout << "usage: printSIFFile <sif file> <bool: print history> <bool: print coeffs>" << std::endl;
  else {
    std::string filename = argv[1];
    int printHistory = atoi(argv[2]);
    int printCoeffs = atoi(argv[3]);

    // open sif file
    SIFFile sfile(filename);
    // always print out infos stored in header
    sfile.printHeader();
    // if required, print out decomposition history
    if (printHistory == 1)
      sfile.printHistory();
    // if required, print out shapelet coeffs (+ errors)
    if (printCoeffs == 1)
      sfile.printCoefficients();
  }
}
