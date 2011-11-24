#ifndef SHAPELENS_IO_H
#define SHAPELENS_IO_H

#include <fitsio.h>
#include <iostream>
#include <string>
#include <map>
#include <stdexcept>
#include <numla/NumMatrix.h>
#include <boost/lexical_cast.hpp>
#include <gsl/gsl_rng.h>
#include "../Typedef.h"
#include "../ShapeLensConfig.h"
#include "../frame/Grid.h"
#include "../frame/WCSTransformation.h"

namespace shapelens {
template <class T>
class Image;

/// Class for FITS-related functions.
/// \b NOTE: The functions will throw a instance of \p std::exception 
/// in case of a failure.
class IO {
  // helper function to get pointer to data section fits functions
  template<class T>
    inline static  T* data(T& val) { 
    return &val;
  }    
  template<class T>
    inline static  T* data(Point<T>& p) {
    return p.c_array();
  }
  inline static char* data(std::string s) {
    return const_cast<char*>(s.c_str());
  }
  
 public:
 
  /// Return FITS Image format definition, based on the type of \p entry.
  /// - <tt>char, unsigned char -> BYTE_IMG</tt>
  /// - <tt>int -> SHORT_IMG</tt>
  /// - <tt>unsigned int -> USHORT_IMG</tt>
  /// - <tt>long -> LONG_IMG</tt>
  /// - <tt>unsigned long -> ULONG_IMG</tt>
  /// - <tt>float -> FLOAT_IMG</tt>
  /// - <tt>double -> DOUBLE_IMG</tt>
  template <typename T> inline static
     int getFITSImageFormat(const T& entry) {
    // default type, uses template specialization for other types
    // see below
    return BYTE_IMG;
  }

  /// Return FITS Image format definition, based on the type of \p entry.
  /// - <tt>bool, char, unsigned char -> TBYTE</tt>
  /// - <tt>int -> TINT</tt>
  /// - <tt>unsigned int -> TUINT</tt>
  /// - <tt>long -> TLONG</tt>
  /// - <tt>unsigned long -> TULONG</tt>
  /// - <tt>float -> TFLOAT</tt>
  /// - <tt>double -> TDOUBLE</tt>
  /// - <tt>std::complex<float> -> TCOMPLEX</tt>
  /// - <tt>std::complex<double> -> TDBLCOMPLEX</tt>
  /// - <tt>std::string -> TSTRING</tt>
  template <typename T> inline static
     int getFITSDataType(const T& entry) {
    // default type, uses template specialization for other types
    // see below
    return TBYTE;
  }

  /// Open FITS file.
  /// If <tt>write == false</tt>, the file will be opened in read-only mode.
  static fitsfile* openFITSFile(const std::string& filename, bool write=false);
  /// Open FITS file and go to first table extension.
  /// If <tt>write == false</tt>, the file will be opened in read-only mode.
  static fitsfile* openFITSTable(const std::string& filename, bool write=false);
  /// Create new FITS file.
  /// If the file \p filename already exists, it will be overwritten.
  static fitsfile* createFITSFile(const std::string& filename);
  /// Close FITS file pointer.
  static void closeFITSFile(fitsfile* fptr);
  /// Move to extension \p i (starting from 1) in FITS file.
  static void moveToFITSExtension(fitsfile* fptr, unsigned int i);
  /// Move to extension \p name in FITS file.
  static void moveToFITSExtension(fitsfile* fptr, const std::string& name);
  /// Append \p history to FITS header histroy.
  static void appendFITSHistory(fitsfile *fptr, const std::string& history);
  /// Read FITS keyword cards directly.
  static void readFITSKeyCards(fitsfile *fptr, const std::string& key, std::string& value);
  /// Get name of FITS file from its pointer.
  static std::string getFITSFileName(fitsfile *fptr);

    /// Read in keyword from FITS header.
  template <class T>
    static void readFITSKeyword(fitsfile *fptr, const std::string& key, T& val) {
    int status = 0;
    char* comment = NULL;
    fits_read_key (fptr,getFITSDataType(val), const_cast<char*>(key.c_str()),data(val),comment, &status);
    // FIXME: for whatever reason, the exception below creates a SEGFAULT!
    /*
    if (status != 0) {
      std::ostringstream note;
      note << "IO: Cannot read FITS keyword " << key << " from " << getFITSFileName(fptr);
      throw std::invalid_argument(note.str());
    }
    */
  }

    /// Set/update keyword in FITS file header.
  template <class T>
    static void updateFITSKeyword(fitsfile *fptr, const std::string& keyword, const T& value_, std::string comment = "") {
    T& value = const_cast<T&>(value_);
    int status = 0;
    fits_write_key (fptr, getFITSDataType(value), const_cast<char *>(keyword.c_str()), data(value), comment.c_str(), &status);
    /*
    if (status != 0) {
      std::ostringstream note;
      note << "IO: Cannot update FITS keyword " << keyword << " = " << value_ << " in " << getFITSFileName(fptr);
      throw std::runtime_error(note.str());
      }
    */
  }


  /// Write FITS image from an Image<T>.
  /// The datatype will be automatically adjusted, based on the
  /// result of getFITSImageFormat() and getFITSDataType().
  /// If \p extname is non-empty, the keyword \p EXTNAME is set to \p extname.
  template <class T>
    static void writeFITSImage(fitsfile *fptr, const Image<T>& image, std::string extname="") {
    int dim0 = image.grid.getSize(0);
    int dim1 = image.grid.getSize(1);
    long naxis = 2;      
    long naxes[2] = { dim0, dim1 };
    long npixels = dim0*dim1;

    // define image format and dataformat according to cfitsio definitions
    int imageformat = getFITSImageFormat(image(0));
    int datatype = getFITSDataType(image(0));
    // create HDU
    int status = 0;
    fits_create_img(fptr, imageformat, naxis, naxes, &status);
    // write pixel data
    long firstpix[2] = {1,1};
    fits_write_pix(fptr,datatype,firstpix,npixels,const_cast<T *>(image.c_array()), &status);
    // insert creator and extname keywords
    if (extname != "")
      updateFITSKeyword(fptr, "EXTNAME", extname);
    updateFITSKeyword(fptr, "CREATOR", std::string("ShapeLens++"));
    if (status != 0)
      throw std::runtime_error("IO: Cannot write FITS image in " + getFITSFileName(fptr));
  }
  
   template <class T>
    static void writeFITSImage(fitsfile *fptr, const Image<std::complex<T> >& image, std::string extname="") {
    int dim0 = image.grid.getSize(0);
    int dim1 = image.grid.getSize(1);
    long naxis = 2;      
    long naxes[2] = { dim0, dim1 };
    long npixels = dim0*dim1;
    long firstpix[2] = {1,1};

    // define image format and dataformat according for each component
    T val;
    int imageformat = getFITSImageFormat(val);
    int datatype = getFITSDataType(val);
    NumVector<T> component(image.size());
    // create HDU
    int status = 0;
    fits_create_img(fptr, imageformat, naxis, naxes, &status);
    // write first component
    for(unsigned long i=0; i < component.size(); i++)
      component(i) = real(image(i));
    fits_write_pix(fptr,datatype,firstpix,npixels,const_cast<T *>(component.c_array()), &status);

    // insert creator and extname keywords
    if (extname != "")
      updateFITSKeyword(fptr, "EXTNAME", extname);
    updateFITSKeyword(fptr, "CREATOR", std::string("ShapeLens++"));

    // write second component
    fits_create_img(fptr, imageformat, naxis, naxes, &status);
    for(unsigned long i=0; i < component.size(); i++)
      component(i) = imag(image(i));
    fits_write_pix(fptr,datatype,firstpix,npixels,const_cast<T *>(component.c_array()), &status);
    if (status != 0)
      throw std::runtime_error("IO: Cannot write FITS image to " + getFITSFileName(fptr));
  }

  /// Write FITS image from an NumMatrix<T>.
  /// The datatype will be automatically adjusted, based on the
  /// result of getFITSImageFormat() and getFITSDataType().
  /// If \p extname is non-empty, the keyword \p EXTNAME is set to \p extname.
  template <class T>
    static void writeFITSImage(fitsfile *fptr, const NumMatrix<T>& M, std::string extname="") {
    int dim0 = M.getColumns();
    int dim1 = M.getRows();
    long naxis = 2;      
    long naxes[2] = { dim0, dim1 };
    long npixels = dim0*dim1;

    // define image format and dataformat according to cfitsio definitions
    int imageformat = getFITSImageFormat(M(0,0));
    int datatype = getFITSDataType(M(0,0));
    // create HDU
    int status = 0;
    fits_create_img(fptr, imageformat, naxis, naxes, &status);
    // write pixel data
    long firstpix[2] = {1,1};
    fits_write_pix(fptr,datatype,firstpix,npixels,const_cast<T *>(M.c_array()), &status);
    // insert creator and extname keywords
    if (extname != "")
      updateFITSKeyword(fptr, "EXTNAME", extname);
    updateFITSKeyword(fptr, "CREATOR", std::string("ShapeLens++"));
    if (status != 0)
      throw std::runtime_error("IO: Cannot write FITS image to " + getFITSFileName(fptr));
  }

  /// Read FITS image into NumMatrix<T>.
  /// \p M is adjusted to hold the contents of the image; the image value are 
  /// automatically casted to the type \p T of \p M.
  template <class T>
    static void readFITSImage(fitsfile *fptr, NumMatrix<T>& M) {
    int naxis, status = 0;
    fits_get_img_dim(fptr, &naxis, &status);
    if (naxis!=2)
      throw std::invalid_argument("IO: naxis != 2. Pointer of " + getFITSFileName(fptr) + " does not provide image");

    long naxes[2] = {1,1};
    fits_get_img_size(fptr, naxis, naxes, &status);
    M.resize(naxes[1],naxes[0]);
    long firstpix[2] = {1,1};
    T val;
    int imageformat = getFITSImageFormat(val);
    int datatype = getFITSDataType(val);
    fits_read_pix(fptr, datatype, firstpix, naxes[0]*naxes[1], NULL, M.c_array(), NULL, &status);
    if (status != 0)
      throw std::runtime_error("IO: Cannot read FITS image in " + getFITSFileName(fptr));
  }

  /// Read FITS image into Image<T>.
  /// \p im is adjusted to hold the contents of the image; the image value are 
  /// automatically casted to the type \p T of \p im.\n
  /// \p Image<T>::grid is set to Grid(0,0,N,M), where \p N and \p M 
  /// are the row and column numbers of the FITS image.
  template <class T>
    static void readFITSImage(fitsfile *fptr, Image<T>& im) {
    int naxis, status = 0;
    fits_get_img_dim(fptr, &naxis, &status);
    if (naxis!=2)
      throw std::invalid_argument("IO: naxis != 2. Pointer of " + getFITSFileName(fptr) + " does not provide image");

    long naxes[2] = {1,1};
    fits_get_img_size(fptr, naxis, naxes, &status);
    im.grid.setSize(0,0,naxes[0],naxes[1]);
    // set wcs from fits file if requested
    if (ShapeLensConfig::USE_WCS) {
#ifdef HAS_WCSLIB
      im.grid.setWCS(WCSTransformation(fptr));
#else
      throw std::runtime_error("IO: WCS usage requested, but HAS_WCSToolsLib not specified");
#endif
    }

    im.resize(im.grid.size());
    long firstpix[2] = {1,1};
    T val;
    int datatype = getFITSDataType(val);
    fits_read_pix(fptr, datatype, firstpix, im.size(), NULL, im.c_array(), NULL, &status);
    if (status != 0)
      throw std::runtime_error("IO: Cannot read FITS image from "+ getFITSFileName(fptr));
  }

  template <class T>
   static void readFITSImage(fitsfile *fptr, Image<std::complex<T> >& im) {
    int naxis, status = 0;
    fits_get_img_dim(fptr, &naxis, &status);
    if (naxis!=2)
      throw std::invalid_argument("IO: naxis != 2. Pointer of " + getFITSFileName(fptr) + " does not provide image");

    long naxes[2] = {1,1};
    fits_get_img_size(fptr, naxis, naxes, &status);
    im.grid.setSize(0,0,naxes[0],naxes[1]);
    // set WCS from fits file if requested
    if (ShapeLensConfig::USE_WCS) {
#ifdef HAS_WCSLIB
      im.grid.setWCS(WCSTransformation(fptr));
#else
      throw std::runtime_error("IO: WCS usage requested, but HAS_WCSToolsLib not specified");
#endif
    }

    im.resize(im.grid.size());
    long firstpix[2] = {1,1};
    T val;
    int datatype = getFITSDataType(val);
    NumVector<T> component(im.size());
    fits_read_pix(fptr, datatype, firstpix, im.size(), NULL, component.c_array(), NULL, &status);
    // copy 1st component
    for(unsigned long i=0; i < component.size(); i++)
      im(i) = component(i);
    int hdutype;
    fits_movrel_hdu(fptr, 1, &hdutype, &status);
    fits_read_pix(fptr, datatype, firstpix, im.size(), NULL, component.c_array(), NULL, &status);
    for(unsigned long i=0; i < component.size(); i++)
      im(i) += std::complex<T>(0,component(i));
    if (status != 0)
      throw std::runtime_error("IO: Cannot read FITS image from " + getFITSFileName(fptr));
  }


  /// Get number of rows in FITS table.
  static long getFITSTableRows(fitsfile* fptr);
  /// Get column number of a FITS table with given \p name.
  /// \p name can contain \p * wildcards and is treated case-insensitive.
  static int getFITSTableColumnNumber(fitsfile* fptr, const std::string& name);
  /// Get data type of column \p colnr from a FITS table.
  static int getFITSTableColumnType(fitsfile* fptr, int colnr);
  /// Read \p val from FITS table at \p row and \p colnr.
  /// If a NULL value was stored at this position, \p nullvalue will be 
  /// returned insted.\n
  /// \p row numbers start with 0 according to regular C-style iterations. 
  template <class T>
  static void readFITSTableValue(fitsfile* fptr, long row, int colnr, T& val, T nullvalue = 0) {
    int status = 0, anynull;
    fits_read_col(fptr, getFITSDataType(val), colnr, row+1, 1, 1, &nullvalue, &val, &anynull, &status);
    if (status != 0) {
      std::ostringstream note;
      note << "IO: Cannot read value in (row/col) = (" << row << "/" << colnr << ") from FITS table in " << getFITSFileName(fptr);
      throw std::runtime_error(note.str());
    }
  }

  /// Write PPM file from data on the given grid.
  /// Colorscheme is one of the following:
  /// - "SPECTRUM"
  /// - "WARM"
  /// - "RED"
  /// - "BLUE"
  /// - "GREEN"
  /// - "GRAY"
  /// Scaling is one of the following:
  /// - "LINEAR"
  /// - "SQUARE_ROOT"
  /// - "LOGAITHMIC"
  /// 
  /// \p min and \p max indicate the ends of the accepted range of values, 
  /// values smaller (larger) than <tt>min (max)</tt> are set to 
  /// <tt>min (max)</tt>.
   void writePPMImage(const std::string& filename, const std::string& colorscheme, const std::string& scaling, data_t min, data_t max, const Grid& grid, const NumVector<data_t>& data);
  /// Create RGB representation of data.
  /// This function is usefull for manipulations of the data in RGB space. 
  /// See writePPMImage() for details.
  static void makeRGBImage(NumMatrix<unsigned int>& rgbImage, const std::string& colorscheme, const std::string& scaling, data_t min, data_t max, const Grid& grid, const NumVector<data_t>& data);
  /// Write RGBImage from makeRGBImage() to PPM file.
  static void writeRGB2PPMImage (const std::string& filename, const Grid& grid, const NumMatrix<unsigned int>& rgbImage);
  /// Add uniform noise within noisemean to noisemean+noiselimit.
  static void addUniformNoise(NumVector<data_t>& data, data_t noisemean, data_t noiselimit);
  /// Add uniform noise from preconfigured RNG within noisemean to noisemean+noiselimit.
  static void addUniformNoise(NumVector<data_t>& data, const gsl_rng* r, data_t noisemean, data_t noiselimit);
  /// Add gaussian noise.
  static void addGaussianNoise(NumVector<data_t>& data, data_t noisemean, data_t noisesigma);
  /// Add gaussian noise from preconfigured RNG.
  static void addGaussianNoise(NumVector<data_t>& data, const gsl_rng* r, data_t noisemean, data_t noisesigma);
  /// Add poissonian noise.
  /// The noise of the \f$i\f$th pixel is drawn from a gaussian distribution with
  /// \f$ \sigma_i = \sqrt{\mu_n + data_i}\f$.
  static void addPoissonianNoise(NumVector<data_t>& data, data_t noisemean);
  /// Add poissonian noise from preconfigured RNG.
  static void addPoissonianNoise(NumVector<data_t>& data, const gsl_rng* r, data_t noisemean);
  /// Convolve input with a 3x3 Gaussian.
  static void convolveGaussian(const NumVector<data_t>& input, NumVector<data_t>& result, int width,int height);
 private:
  static int makeColorMatrix(NumMatrix<unsigned int>& m, const std::string& colorscheme);
  static unsigned int getScaledValue(data_t value, int maxcolors, data_t min, data_t max, char scaling);
};

// template specializations
  template<> inline
    int IO::getFITSImageFormat<unsigned int>(const unsigned int& entry) {
    return USHORT_IMG;
  }
  template<> inline
    int IO::getFITSImageFormat<long>(const long& entry) {
    return LONG_IMG;
  }
  template<> inline
    int IO::getFITSImageFormat<unsigned long>(const unsigned long& entry) {
    return SHORT_IMG;
  }
  template<> inline
    int IO::getFITSImageFormat<float>(const float& entry) {
    return FLOAT_IMG;
  }
  template<> inline
    int IO::getFITSImageFormat<double>(const double& entry) {
    return DOUBLE_IMG;
  }
  

  template<> inline
    int IO::getFITSDataType<int>(const int& entry) {
    return TINT;
  }
  template<> inline
    int IO::getFITSDataType<unsigned int>(const unsigned int& entry) {
    return TUINT;
  }
  template<> inline
    int IO::getFITSDataType<long>(const long& entry) {
    return TLONG;
  }
  template<> inline
    int IO::getFITSDataType<unsigned long>(const unsigned long& entry) {
    return TULONG;
  }
  template<> inline
    int IO::getFITSDataType<float>(const float& entry) {
    return TFLOAT;
  }
  template<> inline
    int IO::getFITSDataType<double>(const double& entry) {
    return TDOUBLE;
  }
  template<> inline
    int IO::getFITSDataType<std::complex<float> >(const std::complex<float>& entry) {
    return TCOMPLEX;
  }
  template<> inline
    int IO::getFITSDataType<std::complex<double> >(const std::complex<double>& entry) {
    return TDBLCOMPLEX;
  }
  template<> inline
    int IO::getFITSDataType<Point<float> >(const Point<float>& entry) {
    return TCOMPLEX;
  }
  template<> inline
    int IO::getFITSDataType<Point<double> >(const Point<double>& entry) {
    return TDBLCOMPLEX;
  }
  template<> inline
    int IO::getFITSDataType<std::string>(const std::string& entry) {
    return TSTRING;
  }

  template <> inline
    void IO::readFITSKeyword<std::string>(fitsfile *fptr, const std::string& key, std::string& val) {
    int status = 0;
    char value[FLEN_CARD];
    fits_read_key (fptr,IO::getFITSDataType(val), const_cast<char *>(key.c_str()),&value, NULL, &status);
    val = std::string(value);
  }

} // end namespace
#endif
