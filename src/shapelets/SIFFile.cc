#include <shapelets/SIFFile.h>

using namespace std;

SIFFile::SIFFile(std::string infilename) {
  filename = infilename;
}

void SIFFile::save(std::string historyString, const NumMatrix<double>& cartesianCoeffs, const NumMatrix<double>& errors, const Grid& grid, double beta, const Point2D& xcentroid, double chi2, char fitsFlag, char decompositionFlag, bool regularize, double R) {
    // version of the sif format definition
  int version = 0;
  const char* historyChar = historyString.c_str();
  int lines = cartesianCoeffs.getRows();
  int columns = cartesianCoeffs.getColumns();
  int errorLines = errors.getRows();
  int errorColumns = errors.getColumns();
  int dimension = lines*columns;
  int errorDimension = errorLines*errorColumns;
  const char* file = filename.c_str();
  sifHeader header = { version, lines, columns, errorLines, errorColumns, beta, xcentroid(0), xcentroid(1), grid.getStartPosition(0), grid.getStopPosition(0), grid.getStepsize(0), grid.getStartPosition(1), grid.getStopPosition(1), grid.getStepsize(1), chi2, (int) fitsFlag, (int) decompositionFlag, regularize, R, historyString.length()};
  
  fstream binary_file(file,ios::out|ios::binary);
  saveHeader(binary_file,header);

  double dataArray[dimension], errorArray[errorDimension];
  for (int i = 0; i < lines; i++)
    for (int j = 0; j < columns; j++)
      dataArray[i*columns + j] = cartesianCoeffs(i,j);
  for (int i = 0; i < errorLines; i++)
    for (int j = 0; j < errorColumns; j++)
      errorArray[i*errorColumns + j] = errors(i,j);

  binary_file.write(historyChar,historyString.length()*sizeof(char));
  binary_file.write(reinterpret_cast<char *>(&dataArray),dimension*sizeof(double));
  binary_file.write(reinterpret_cast<char *>(&errorArray),errorDimension*sizeof(double));
  binary_file.close();
}

void SIFFile::saveHeader (fstream& binary_file, sifHeader& header) {
  binary_file.write(reinterpret_cast<char *>(&header.version),sizeof(int));
  binary_file.write(reinterpret_cast<char *>(&header.lines),sizeof(int));
  binary_file.write(reinterpret_cast<char *>(&header.columns),sizeof(int));
  binary_file.write(reinterpret_cast<char *>(&header.errorLines),sizeof(int));
  binary_file.write(reinterpret_cast<char *>(&header.errorColumns),sizeof(int));
  binary_file.write(reinterpret_cast<char *>(&header.beta),sizeof(double));
  binary_file.write(reinterpret_cast<char *>(&header.xcentroid0),sizeof(double));
  binary_file.write(reinterpret_cast<char *>(&header.xcentroid1),sizeof(double));
  binary_file.write(reinterpret_cast<char *>(&header.gridstart0),sizeof(double));
  binary_file.write(reinterpret_cast<char *>(&header.gridstop0),sizeof(double));
  binary_file.write(reinterpret_cast<char *>(&header.gridstepsize0),sizeof(double));
  binary_file.write(reinterpret_cast<char *>(&header.gridstart1),sizeof(double));
  binary_file.write(reinterpret_cast<char *>(&header.gridstop1),sizeof(double));
  binary_file.write(reinterpret_cast<char *>(&header.gridstepsize1),sizeof(double));
  binary_file.write(reinterpret_cast<char *>(&header.chi2),sizeof(double));
  binary_file.write(reinterpret_cast<char *>(&header.fitsFlag),sizeof(int));
  binary_file.write(reinterpret_cast<char *>(&header.decompositionFlag),sizeof(int));
  binary_file.write(reinterpret_cast<char *>(&header.regularize),sizeof(int));
  binary_file.write(reinterpret_cast<char *>(&header.R),sizeof(double));
  binary_file.write(reinterpret_cast<char *>(&header.historylength),sizeof(int));
}
  
void SIFFile::load(std::string& historyString, NumMatrix<double>& cartesianCoeffs, NumMatrix<double>& errors, Grid& grid, double& beta, Point2D& xcentroid, double& chi2, char& fitsFlag, char& decompositionFlag, bool& regularize, double& R) {
  sifHeader header;

  const char* file = filename.c_str();
  fstream binary_file (file,ios::binary|ios::in);
  // check if file handle is valid
  testFileHandle(binary_file);
  binary_file.clear();

  // first read in header
  loadHeader(binary_file,header);

  // we have to expand historyChar by one 0 char, to finalize the string
  char historyChar[header.historylength+1];
  binary_file.read(historyChar,header.historylength*sizeof(char));
  historyChar[header.historylength]=0;
  historyString = string(historyChar);
  fitsFlag = (char) header.fitsFlag;
  decompositionFlag = (char) header.decompositionFlag;
  regularize = (bool) header.regularize;
  R = header.R;

  int dimension = header.lines * header.columns;
  int errorDimension = header.errorLines * header.errorColumns;
  double dataArray[dimension], errorArray[errorDimension];
  binary_file.read(reinterpret_cast<char *>(&dataArray),dimension*sizeof(double));
  binary_file.read(reinterpret_cast<char *>(&errorArray),errorDimension*sizeof(double));
  binary_file.close();

  cartesianCoeffs = NumMatrix<double>(header.lines, header.columns);
  errors = NumMatrix<double>(header.errorLines, header.errorColumns);
  for (int i = 0; i < header.lines; i++)
    for (int j = 0; j < header.columns; j++)
      cartesianCoeffs(i,j) = dataArray[i*header.columns + j];
  for (int i = 0; i < header.errorLines; i++)
    for (int j = 0; j < header.errorColumns; j++)
      errors(i,j) = errorArray[i*header.errorColumns + j];

  
  grid = Grid(header.gridstart0, header.gridstop0, header.gridstepsize0, header.gridstart1, header.gridstop1, header.gridstepsize1);
  beta = header.beta;
  xcentroid = Point2D(header.xcentroid0,header.xcentroid1);
  chi2 = header.chi2;
}

void SIFFile::loadHeader(sifHeader& header) {
  const char* file = filename.c_str();
  fstream binary_file (file,ios::binary|ios::in);
  testFileHandle(binary_file);
  binary_file.clear();
  loadHeader(binary_file,header);
  binary_file.close();
}

void SIFFile::loadHeader(fstream& binary_file, sifHeader& header) {
  binary_file.read(reinterpret_cast<char *>(&header.version),sizeof(int));
  binary_file.read(reinterpret_cast<char *>(&header.lines),sizeof(int));
  binary_file.read(reinterpret_cast<char *>(&header.columns),sizeof(int));
  binary_file.read(reinterpret_cast<char *>(&header.errorLines),sizeof(int));
  binary_file.read(reinterpret_cast<char *>(&header.errorColumns),sizeof(int));
  binary_file.read(reinterpret_cast<char *>(&header.beta),sizeof(double));
  binary_file.read(reinterpret_cast<char *>(&header.xcentroid0),sizeof(double));
  binary_file.read(reinterpret_cast<char *>(&header.xcentroid1),sizeof(double));
  binary_file.read(reinterpret_cast<char *>(&header.gridstart0),sizeof(double));
  binary_file.read(reinterpret_cast<char *>(&header.gridstop0),sizeof(double));
  binary_file.read(reinterpret_cast<char *>(&header.gridstepsize0),sizeof(double));
  binary_file.read(reinterpret_cast<char *>(&header.gridstart1),sizeof(double));
  binary_file.read(reinterpret_cast<char *>(&header.gridstop1),sizeof(double));
  binary_file.read(reinterpret_cast<char *>(&header.gridstepsize1),sizeof(double));
  binary_file.read(reinterpret_cast<char *>(&header.chi2),sizeof(double));
  binary_file.read(reinterpret_cast<char *>(&header.fitsFlag),sizeof(int));
  binary_file.read(reinterpret_cast<char *>(&header.decompositionFlag),sizeof(int));
  binary_file.read(reinterpret_cast<char *>(&header.regularize),sizeof(int));
  binary_file.read(reinterpret_cast<char *>(&header.R),sizeof(double));
  binary_file.read(reinterpret_cast<char *>(&header.historylength),sizeof(int));
}

void SIFFile::printHeader() {
  sifHeader header;
  const char* file = filename.c_str();
  fstream binary_file (file,ios::binary|ios::in);
  testFileHandle(binary_file);
  binary_file.clear();
  loadHeader(binary_file, header);
  binary_file.close();

  std::cout << "SIF header information" << std::endl;
  std::cout << "Filename:\t" << filename << std::endl;
  std::cout << "SIF version:\t" <<  header.version << std::endl;
  std::cout << "Coefficients:\t" << header.lines <<"x" << header.columns << std::endl;
  std::cout << "Coff. errors:\t" << header.errorLines <<"x" << header.errorColumns << std::endl;;
  std::cout << "Beta:\t\t" << header.beta << std::endl;
  std::cout << "Centroid:\t" << header.xcentroid0 << "/" << header.xcentroid1 << std::endl;
  std::cout << "Grid:\t\t" << header.gridstart0 << ".." << header.gridstop0 << " (" <<  header.gridstepsize0 << "); " << header.gridstart1 << ".." << header.gridstop1 << " (" <<  header.gridstepsize1 << ")" << std::endl;
  std::cout << "Chi^2:\t\t" << header.chi2 << std::endl;
  std::cout << "Flags:\t\t" << header.fitsFlag << "/" << header.decompositionFlag << std::endl;
  std::cout << "Regularized:\t" << header.regularize;
  if (header.regularize)
    std::cout << " (R = " << header.R << ")" << std::endl;
  else 
    std::cout << std::endl;
  std::cout << "History length:\t" << header.historylength << std::endl;
}

void SIFFile::printHistory() {
  std::string historyString;
  NumMatrix<double> cartesianCoeffs,errors;
  Grid grid;
  Point2D xcentroid;
  double beta, chi2, R;
  char fitsFlag, decompositionFlag;
  bool regularized;
  load(historyString,cartesianCoeffs,errors,grid,beta,xcentroid,chi2,fitsFlag,decompositionFlag,regularized,R);

  std::cout << historyString;
}

void SIFFile::printCoefficients() {
  std::string historyString;
  NumMatrix<double> cartesianCoeffs,errors;
  Grid grid;
  Point2D xcentroid;
  double beta, chi2, R;
  char fitsFlag, decompositionFlag;
  bool regularized;
  load(historyString,cartesianCoeffs,errors,grid,beta,xcentroid,chi2,fitsFlag,decompositionFlag,regularized,R);

  std::cout << "Cartesian Coefficients:"<< std::endl;
  std::cout << cartesianCoeffs << std::endl;
  if (errors.getRows() != 0) {
    std::cout << "Cartesian Errors:"<< std::endl;
    std::cout << errors << std::endl;
  }
}

void SIFFile::testFileHandle(fstream& binary_file) {
  if (binary_file.fail()) {
    std::cout << "SIFFile: file does not exist!" << std::endl;
    terminate();
  }
}
