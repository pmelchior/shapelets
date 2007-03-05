#ifndef FITSIMAGE_H
#define FITSIMAGE_H

#include <fitsio.h>
#include <string>
#include <NumVector.h>
#include <Grid.h>
#include <IO.h>

/// Class for manipulation FITS files.
/// The class allows reading the pHDU and extHDUs of standard conform FITS files.
template <class T>
class FitsImage {
 public:
  /// Argumented constructor.
  /// The extension has to given in the standard cfitsio way as
  /// filename[extension].
  FitsImage(std::string infilename) {
    filename = infilename;
    read();
  }
  /// Get axis size of the whole image in given direction.
  unsigned int getSize(bool direction) {
    if (direction == 0) return axsize0;
    else return axsize1;
  }
  /// Get contents of image as NumVector.
  const NumVector<T>& getData() {
    return data;
  }
  /// Access contents of image.
  NumVector<T>& accessData() {
    return data;
  }
  /// Get the Grid for the image.
  /// The grid is defined to range from 0 to getAxisSize(i) in steps of 1.
  const Grid& getGrid() {
    return grid;
  }
  /// Return the number of pixels in the image.
  unsigned int getNumberOfPixels() {
    return axsize0*axsize1;
  }
  /// Set the x and y coordinate from the pixel number.
  /// The coordinate system is defined such, that the left lower corder of the image
  /// has the coordinates (0,0) and each pixel has unit length.
  void getCoords(unsigned int pixel, int& x, int& y) {
    x = pixel%axsize0;
    y = pixel/axsize0;
  }
  /// Get the pixel number from the x and y corrdinates.
  unsigned int getPixel(int x, int y) {
    return (unsigned int) x + y*axsize0;
  }
  /// Set the filename.
  /// When this method is called with a different filename, 
  /// the specified file will be opened and read.
  void setFilename(std::string infilename) {
    filename = infilename;
    read();
  }
  /// Get the filename of the Fits file.
  std::string getFilename() {
    return filename;
  }

 private:
  unsigned int axsize0, axsize1;
  NumVector<T> data;
  Grid grid;
  std::string filename;
  void read() {
    fitsfile *fptr;
    int status = 0;
    fits_open_file(&fptr, filename.c_str(), READONLY, &status);
    int naxis;
    fits_get_img_dim(fptr, &naxis, &status);
    if (naxis!=2) {
      std::cout << "FitsImage: naxis != 2. This is not a FITS image!" << std::endl;
      std::terminate();
    } else {
      long naxes[2] = {1,1};
      fits_get_img_size(fptr, naxis, naxes, &status);
      axsize0 = naxes[0];
      axsize1 = naxes[1];
      grid = Grid(0,axsize0-1,1,0,axsize1-1,1);
      long npixels = axsize0*axsize1;
      data.resize(npixels);
      long firstpix[2] = {1,1};
      T val;
      int imageformat, datatype;
      setFITSTypes(val,imageformat,datatype);
      fits_read_pix(fptr, datatype, firstpix, npixels, NULL,data.c_array(), NULL, &status);
      fits_close_file(fptr, &status);
    }
  }
};

#endif
