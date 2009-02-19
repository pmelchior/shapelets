#ifndef IMAGE_H
#define IMAGE_H

#include <fitsio.h>
#include <string>
#include <NumVector.h>
#include <Typedef.h>
#include <History.h>
#include <IO.h>
#include <frame/Grid.h>

/// Class for image data.
/// Internally, the class represents image data as a vector. But the interface
/// allows two-dimensional access.\n
/// The class can read in data directly from and save it to FITS files.

template <class T>
class Image : public NumVector<T> {
  typedef NumVector<T> data;
 public:
  /// Default constructor.
  Image() : NumVector<T>() {
    filename = "";
    history.clear();
  }
  /// Constructor with given image size \f$N\times M\f$.
  Image(unsigned int N, unsigned int M) : NumVector<T>(N*M) {
    filename = "";
    history.clear();
    grid = Grid(0,0,N,M);
  }
  /// Argumented constructor for reading from a FITS file.
  /// The extension can be given in the standard cfitsio way as
  /// filename[extension].
  Image(std::string infilename) : NumVector<T>() {
    filename = infilename;
    history.clear();
    read();
  }
  /// Copy operator from base-class.
  void operator=(const NumVector<T>& v) {
    NumVector<T>::operator=(v);
  }
  /// Access operator using pixel index.
  T& operator()(unsigned long i) {
    return NumVector<T>::operator()(i);
  }
  const T& operator()(unsigned long i) const {
    return NumVector<T>::operator()(i);
  }
  /// Access operator using pixel coordinates.
  T& operator()(unsigned long x, unsigned long y) {
    return NumVector<T>::operator()(y*Image::getSize(0) + x);
  }
  /// const access operator using pixel coordinates.
  const T& operator()(unsigned long x, unsigned long y) const {
    return NumVector<T>::operator()(y*Image::getSize(0) + x);
  }
  NumVector<T>& accessNumVector() {
    return *this;
  }
  const NumVector<T>& getNumVector() const {
    return *this;
  }
  /// Slice a sub-image, specified by edge-points \p P1 and \p P2.
  Image<T> slice(const Point2D<grid_t>& P1, const Point2D<grid_t>& P2) const {
    int xmin = P1(0), xmax = P2(0), ymin = P1(1), ymax = P2(1);
    int axis0 = xmax-xmin;
    Image<T> sub;
    sub.resize((xmax-xmin)*(ymax-ymin));

    // lop over all object pixels
    for (int i =0; i < sub.size(); i++) {
      // old coordinates derived from new pixel index i
      int x = i%axis0 + xmin;
      int y = i/axis0 + ymin;
      if (x>=0 && x < Image::getSize(0) && y >= 0 && y < Image::getSize(1))
	sub(i) = Image<T>::operator()(x,y);
    }
    sub.grid = Grid(xmin,ymin,xmax-xmin,ymax-ymin);
    sub.filename = filename;
    sub.history.setSilent();
    sub.history << "# Slice from " << filename << " in the area (";
    sub.history << xmin << "/" << ymin << ") -> (";
    sub.history << xmax << "/" << ymax << ")" << std::endl;
    return sub;
  }
  /// Get axis size of the whole image in given direction.
  unsigned int getSize(bool direction) const {
    return grid.getSize(direction);
  }
  /// Get the filename of the Fits file.
  std::string getFilename() const {
    return filename;
  }
  /// Save as FITS file.
  void save(std::string filename) const {
    fitsfile* fptr = IO::createFITSFile(filename);
    IO::writeFITSImage(fptr,grid,*this);
    // if history is not empty, append history to FITS header
    if (!history.isEmpty())
      IO::appendFITSHistory(fptr,history.str());
    IO::closeFITSFile(fptr);
  }
  /// The Grid this Image is defined on.
  /// The grid is defined to range from 0 to getAxisSize(i)-1 in steps of 1.
  Grid grid;
  /// The image history.
  History history;
  
  // Legacy function
  const Grid& getGrid() const {
    return grid;
  }

 protected:
  std::string filename;
 private:
  void read() {
    fitsfile *fptr = IO::openFITSFile(filename);
    int status = IO::readFITSImage(fptr,grid,*this);
    status = IO::closeFITSFile(fptr);
  }
};
#endif
