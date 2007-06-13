#ifndef IMAGE_H
#define IMAGE_H

#include <fitsio.h>
#include <string>
#include <NumVector.h>
#include <Grid.h>
#include <IO.h>

/// Class for image data.
/// The class represents image data as a vector.\n
/// It is particularly suited for reading in data from FITS files.

template <class T>
class Image : public NumVector<T> {
  typedef NumVector<T> data;
 public:
  /// Default constructor.
  /// Only usefull when constructing a FitsImage from data in memory.
  Image() : NumVector<T>() {
    filename = "";
  }
  /// Argumented constructor for reading from a FITS file.
  /// The extension has to given in the standard cfitsio way as
  /// filename[extension].
  Image(std::string infilename) : NumVector<T>() {
    filename = infilename;
    read();
  }
  /// Get axis size of the whole image in given direction.
  unsigned int getSize(bool direction) const {
    if (direction == 0) return grid.getSize(0);
    else return grid.getSize(1);
  }
  /// Get contents of image as NumVector.
  const NumVector<T>& getData() const {
    //return data;
    return *this;
  }
  /// Access contents of image.
  NumVector<T>& accessData() {
    //    return data;
    return *this;
  }
  /// Get the Grid for the image.
  /// The grid is defined to range from 0 to getAxisSize(i)-1 in steps of 1.
  const Grid& getGrid() const {
    return grid;
  }
  /// Access the Grid of the image.
  Grid& accessGrid() {
    return grid;
  }
  /// Get the filename of the Fits file.
  std::string getFilename() const {
    return filename;
  }
  /// Save as FITS file.
  void save(std::string fitsfile) const {
    writeFITSFile(fitsfile,grid,*this);
  }

 private:
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
      unsigned int axsize0, axsize1;
      axsize0 = naxes[0];
      axsize1 = naxes[1];
      grid = Grid(0,axsize0-1,1,0,axsize1-1,1);
      long npixels = axsize0*axsize1;
      data::resize(npixels);
      long firstpix[2] = {1,1};
      T val;
      int imageformat, datatype;
      setFITSTypes(val,imageformat,datatype);
      fits_read_pix(fptr, datatype, firstpix, npixels, NULL,data::c_array(), NULL, &status);
      fits_close_file(fptr, &status);
    }
  }
};

#endif
