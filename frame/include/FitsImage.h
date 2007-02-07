#ifndef FITSIMAGE_H
#define FITSIMAGE_H

#include <CCfits/CCfits>
#include <string>
#include <NumVector.h>
#include <Grid.h>

/// Class for manipulation FITS files.
/// The class allows reading the pHDU and extHDUs of standard conform FITS files.
/// \todo Reading as template

class FitsImage {
 public:
  /// Default constructor.
  FitsImage();
  /// Argumented constructor.
  FitsImage(string filename);
  /// Read Fits image filename, using given extension.
  /// Extension = 0 uses pHDU, Extension > 0 uses extHDU.
  void read(int ext);
  /// Print image header to stdout.
  void printHeader();
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
  /// For use with the default constructor only.
  void setFilename(std::string filename);
  /// The filename of the Fits file.
  std::string filename;

 private:
  std::auto_ptr<CCfits::FITS> pInfile;
  unsigned int axsize0, axsize1;
  NumVector<double> data;
  Grid grid;
  template <class C> void getContent (C& image);
  bool isRead;
};

#endif
