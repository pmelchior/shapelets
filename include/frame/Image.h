#ifndef IMAGE_H
#define IMAGE_H

#include <fitsio.h>
#include <string>
#include <NumVector.h>
#include <Typedef.h>
#include <IO.h>
#include <History.h>
#include <frame/Grid.h>

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
    history.clear();
  }
  /// Argumented constructor for reading from a FITS file.
  /// The extension has to given in the standard cfitsio way as
  /// filename[extension].
  Image(std::string infilename) : NumVector<T>() {
    filename = infilename;
    history.clear();
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
    // if history is not empty, append history to FITS header
    if (!history.isEmpty())
      appendFITSHistory(fitsfile,history.str());
  }
  ///  Get image history.
  const History& getHistory() const {
    return history;
  }
  /// Access image history.
  History& accessHistory() {
    return history;
  }

 private:
  Grid grid;
  std::string filename;
  History history;
  void read() {
    fitsfile *fptr = openFITSFile(filename);
    int status = readFITSImage(fptr,grid,*this);
    fits_close_file(fptr, &status);
  }
};

#endif
