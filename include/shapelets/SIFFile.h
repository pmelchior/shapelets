#ifndef SIFFILE_H
#define SIFFILE_H

#include <fstream>
#include <string>
#include <iostream>
#include <NumMatrix.h>
#include <Typedef.h>
#include <frame/Grid.h>
#include <frame/Point2D.h>

/// Define the header of the sif files.
struct sifHeader {
  /// The SIF version: default is 1
  int version;
  /// The number of lines in the cartesian coefficient matrix
  int lines;
  /// The number of columns in the cartesian coefficient matrix
  int columns;
  /// The number of lines in the error matrix.
  /// This is either identical to lines or 0.
  int errorLines;
  /// The number of columns in the error matrix.
  /// This is either identical to columns or 0.
  int errorColumns;
  /// The shapelet scale size
  data_t beta;
  /// 1st component of the centroid position
  data_t xcentroid0;
  /// 2nd component of the centroid position
  data_t xcentroid1;
  /// Starting point of the grid in the first direction
  data_t gridstart0;
  /// Stopping point of the grid in the first direction
  data_t gridstop0;
  /// Stepsize of the grid in the first direction
  data_t gridstepsize0;
  /// Starting point of the grid in the second direction
  data_t gridstart1;
  /// Stopping point of the grid in the second direction 
  data_t gridstop1;
  /// Stepsize of the grid in the first direction
  data_t gridstepsize1;
  /// Chi^2 of the decomposition (or 0 if a artificial image)
  data_t chi2;
  /// The detection flag from Fitsimage
  int fitsFlag;
  /// The decomposition flag from OptimaDecomposition2D
  int decompositionFlag;
  /// Whether Regularization was employed during decomposition
  int regularize;
  /// The \f$R\f$ level reached after regularization
  data_t R;
  /// The length of the image history (in chars)
  int historylength;
};

/// Definition and usage the ShapeletImageFormat.
/// The SIF (Shapelet Image Format) is used to store all information that specifies
/// a ShapeletObject.

class SIFFile {
 public:
  /// Argumented constructor.
  SIFFile(std::string filename);

  /// Save the given information to filename
  void save(std::string historyString, const NumMatrix<data_t>& cartesianCoeffs, const NumMatrix<data_t>& errors, const Grid& grid, data_t beta, const Point2D& xcentroid, data_t chi2, char fitsFlag, char decompositionFlag, bool regularize, data_t R);
  /// Load the shapelet information from filename
  void load(std::string& historyString, NumMatrix<data_t>& cartesianCoeffs, NumMatrix<data_t>& errors, Grid& grid, data_t& beta, Point2D& xcentroid, data_t& chi2, char& fitsFlag, char& decompositionFlag, bool& regularize, data_t& R);
  /// Load the hader information into the sifHeader struct.
  void loadHeader(sifHeader& header);
  /// Print the information in the SIF header to stdout.
  void printHeader();
  /// Print image history to stdout.
  void printHistory();
  /// Print cartesian shapelet coefficients and errors to stdout.
  void printCoefficients();
 private:
  std::string filename;
  void loadHeader(std::fstream& binary_file, sifHeader& header);
  void saveHeader(std::fstream& binary_file, sifHeader& header);
  void testFileHandle(std::fstream& binary_file);
};

#endif
