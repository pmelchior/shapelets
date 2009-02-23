#ifndef IMAGE_H
#define IMAGE_H

#include <fitsio.h>
#include <string>
#include <NumVector.h>
#include <Typedef.h>
#include <utils/History.h>
#include <utils/IO.h>
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
    basefilename = "???";
    history.clear();
  }
  /// Constructor with given image size \f$N\times M\f$.
  Image(unsigned int N, unsigned int M) : NumVector<T>(N*M) {
    basefilename = "???";
    history.clear();
    grid = Grid(0,0,N,M);
  }
  /// Argumented constructor for reading from a FITS file.
  /// The extension can be given in the standard cfitsio way as
  /// \p filename[extension].
  Image(std::string filename) : NumVector<T>() {
    basefilename = filename;
    history.clear();
    read();
  }
  /// Copy constructor from base class.
  Image(const NumVector<T>& v) :  NumVector<T>(v) {
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
  void slice(Image<T>& sub, const Point2D<grid_t>& P1, const Point2D<grid_t>& P2) const {
    int xmin = P1(0), xmax = P2(0), ymin = P1(1), ymax = P2(1);
    int axis0 = xmax-xmin;
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
    sub.basefilename = basefilename;
    sub.history.setSilent();
    sub.history << "# Slice from " << basefilename << " in the area (";
    sub.history << xmin << "/" << ymin << ") -> (";
    sub.history << xmax << "/" << ymax << ")" << std::endl;
  }
  /// Create sampled version if this image by averaging over
  /// \p sampling pixels in each direction
  void sample(Image<T>& sampled, int sampling) {
    T av;
    int highres_x, highres_y;
    if (sampled.size() != Image<T>::size()/(sampling*sampling)) {
      sampled.resize(Image<T>::size()/(sampling*sampling));
      sampled.grid = Grid(0,0,Image<T>::getSize(0)/sampling,Image<T>::getSize(1)/sampling);
    }
    for (int x=0; x < sampled.getSize(0); x++) {
      for (int y=0; y < sampled.getSize(1); y++) {
	av = 0;
	// left-lower starting point for sampling
	if (sampling%2==0) {
	  for (int nx=0; nx<sampling; nx++)
	    for (int ny=0; ny<sampling; ny++)
	      av += Image<T>::operator()(x*sampling + nx, y*sampling + ny);
	} else { // centered sampling
	  for (int nx=-sampling/2; nx<sampling/2; nx++) {
	    for (int ny=-sampling/2; ny<sampling/2; ny++) {
	      highres_x = x*sampling + nx;
	      highres_y = y*sampling + ny;
	      // if pixel get out of bounds, mirror at the center
	      if (highres_x < 0 || highres_x >= Image<T>::getSize(0))
		highres_x = Image<T>::getSize(0) - abs(highres_x);
	      if (highres_y < 0 || highres_y >= Image<T>::getSize(1))
		highres_y = Image<T>::getSize(1) - abs(highres_y);
	      av += Image<T>::operator()(highres_x,highres_y);
	    }
	  }
	}
	sampled(x,y) = av;
      }
    }
    sampled.basefilename = basefilename;
    sampled.history.setSilent();
    sampled.history << "# Downsampled from " << basefilename << " by factor " << sampling << std::endl;
  }
  /// Get axis size of the whole image in given direction.
  unsigned int getSize(bool direction) const {
    return grid.getSize(direction);
  }
  /// Get the filename of the Fits file.
  std::string getFilename() const {
    return basefilename;
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
  /// The file this data is derived from.
  /// The string should not exceed 80 characters, otherwise it becomes 
  /// eventually truncated during save().
  std::string basefilename;
  
  // Legacy function
  const Grid& getGrid() const {
    return grid;
  }
 private:
  void read() {
    fitsfile *fptr = IO::openFITSFile(basefilename);
    int status = IO::readFITSImage(fptr,grid,*this);
    status = IO::closeFITSFile(fptr);
  }
};
#endif
