#ifndef FITSIMAGE_H
#define FITSIMAGE_H

#include <fitsio.h>
#include <string>
#include <NumVector.h>
#include <Grid.h>

/// Class for manipulation FITS files.
/// The class allows reading the pHDU and extHDUs of standard conform FITS files.
/// \todo Reading as template

class FitsImage {
 public:
  /// Argumented constructor.
  /// The extension has to given in the standard cfitsio way as
  /// filename[extension].
  FitsImage(std::string filename);
  /// Get axis size of the whole image in given direction.
  unsigned int getSize(bool direction);
  /// Get contents of image as NumVector.
  const NumVector<double>& getData();
  /// Access contents of image.
  NumVector<double>& accessData();
  /// Get the Grid for the image.
  /// The grid is defined to range from 0 to getAxisSize(i) in steps of 1.
  const Grid& getGrid();
  /// Return the number of pixels in the image.
  unsigned int getNumberOfPixels();
  /// Set the x and y coordinate from the pixel number.
  /// The coordinate system is defined such, that the left lower corder of the image
  /// has the coordinates (0,0) and each pixel has unit length.
  void getCoords(unsigned int pixel, int& x, int& y);
  /// Get the pixel number from the x and y corrdinates.
  unsigned int getPixel(int x, int y);
  /// Set the filename.
  /// When this method is called with a different filename, 
  /// the specified file will be opened and read.
  void setFilename(std::string filename);
  /// Get the filename of the Fits file.
  std::string getFilename();

 private:
  unsigned int axsize0, axsize1;
  NumVector<double> data;
  Grid grid;
  std::string filename;
  void read();
};

#endif
