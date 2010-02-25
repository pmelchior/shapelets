#ifndef SHAPELENS_IO_H
#define SHAPELENS_IO_H

#include <fitsio.h>
#include <iostream>
#include <string>
#include <map>
#include <stdexcept>
#include <numla/NumMatrix.h>
#include "../Typedef.h"
#include "../frame/Grid.h"

namespace shapelens {
template <class T>
class Image;

/// Functions for reading and writing into several formats\n\n
/// \b NOTE: The functions will throw a instance of \p std::exception 
/// in case of a failure.
class IO {
  public:

  /// Return FITS Image format definition, based on the type of \p entry.
  /// - <tt>char, unsigned char -> BYTE_IMG</tt>
  /// - <tt>int -> SHORT_IMG</tt>
  /// - <tt>unsigned int -> USHORT_IMG</tt>
  /// - <tt>long -> LONG_IMG</tt>
  /// - <tt>unsigned long -> ULONG_IMG</tt>
  /// - <tt>float -> FLOAT_IMG</tt>
  /// - <tt>double -> DOUBLE_IMG</tt>
  template <typename T> inline
    static int getFITSImageFormat(const T& entry) {
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
  /// - <tt>complex<float> -> TCOMPLEX</tt>
  /// - <tt>complex<double> -> TDBLCOMPLEX</tt>
  /// - <tt>std::string -> TSTRING</tt>
  template <typename T> inline
    static int getFITSDataType(const T& entry) {
    // default type, uses template specialization for other types
    // see below
    return TBYTE;
  }

  /// Open FITS file.
  /// If <tt>write == false</tt>, the file will be opened in read-only mode.
  static fitsfile* openFITSFile(std::string filename, bool write=false);
  /// Create new FITS file.
  /// If the file \p filename already exists, it will be overwritten.
  static fitsfile* createFITSFile(std::string filename);
  /// Close FITS file pointer.
  static void closeFITSFile(fitsfile* fptr);
  /// Set/update std::string keyword in FITS file header.
  static void updateFITSKeywordString(fitsfile *outfptr, std::string keyword, std::string value, std::string comment="");
  /// Append \p history to FITS header histroy.
  static void appendFITSHistory(fitsfile *outfptr, std::string history);
  /// Read std::string keyword from FITS header.
  static void readFITSKeywordString(fitsfile *fptr, std::string key, std::string& val);
  /// Read FITS keyword cards directly.
  static void readFITSKeyCards(fitsfile *fptr, std::string key, std::string& value);

  /// Write FITS image from an Image<T>.
  /// The datatype will be automatically adjusted, based on the
  /// result of getFITSImageFormat() and getFITSDataType().
  /// If \p extname is non-empty, the keyword \p EXTNAME is set to \p extname.
  template <class T>
    static void writeFITSImage(fitsfile *outfptr, const Image<T>& image, std::string extname="") {
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
    fits_create_img(outfptr, imageformat, naxis, naxes, &status);
    // write pixel data
    long firstpix[2] = {1,1};
    fits_write_pix(outfptr,datatype,firstpix,npixels,const_cast<T *>(image.c_array()), &status);
    // insert creator and extname keywords
    if (extname != "")
      updateFITSKeywordString (outfptr, "EXTNAME", extname);
    updateFITSKeywordString (outfptr, "CREATOR", "ShapeLens++");
    if (status != 0)
      throw std::runtime_error("IO: Cannot write FITS image!");
  }
  
   template <class T>
    static void writeFITSImage(fitsfile *outfptr, const Image<complex<T> >& image, std::string extname="") {
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
    fits_create_img(outfptr, imageformat, naxis, naxes, &status);
    // write first component
    for(unsigned long i=0; i < component.size(); i++)
      component(i) = real(image(i));
    fits_write_pix(outfptr,datatype,firstpix,npixels,const_cast<T *>(component.c_array()), &status);

    // insert creator and extname keywords
    if (extname != "")
      updateFITSKeywordString (outfptr, "EXTNAME", extname);
    updateFITSKeywordString (outfptr, "CREATOR", "ShapeLens++");

    // write second component
    fits_create_img(outfptr, imageformat, naxis, naxes, &status);
    for(unsigned long i=0; i < component.size(); i++)
      component(i) = imag(image(i));
    fits_write_pix(outfptr,datatype,firstpix,npixels,const_cast<T *>(component.c_array()), &status);
    if (status != 0)
      throw std::runtime_error("IO: Cannot write FITS image!");
  }

  /// Write FITS image from an NumMatrix<T>.
  /// The datatype will be automatically adjusted, based on the
  /// result of getFITSImageFormat() and getFITSDataType().
  /// If \p extname is non-empty, the keyword \p EXTNAME is set to \p extname.
  template <class T>
    static void writeFITSImage(fitsfile *outfptr, const NumMatrix<T>& M, std::string extname="") {
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
    fits_create_img(outfptr, imageformat, naxis, naxes, &status);
    // write pixel data
    long firstpix[2] = {1,1};
    fits_write_pix(outfptr,datatype,firstpix,npixels,const_cast<T *>(M.c_array()), &status);
    // insert creator and extname keywords
    if (extname != "")
      updateFITSKeywordString (outfptr, "EXTNAME", extname);
    updateFITSKeywordString (outfptr, "CREATOR", "ShapeLens++");
    if (status != 0)
      throw std::runtime_error("IO: FITS image could not be writted!");
  }

  /// Set/update keyword in FITS file header.
  /// For setting string keywords, use updateFITSKeywordString() instead.
  template <class T>
    static void updateFITSKeyword(fitsfile *outfptr, std::string keyword, T value, std::string comment="") {
    int status = 0;
    fits_write_key (outfptr, getFITSDataType(value), const_cast<char *>(keyword.c_str()), &value, const_cast<char *>(comment.c_str()), &status);
    if (status != 0)
      throw std::runtime_error("IO: Cannot update FITS keyword!");
  }

  /// Read FITS image into NumMatrix<T>.
  /// \p M is adjusted to hold the contents of the image; the image value are 
  /// automatically casted to the type \p T of \p M.
  template <class T>
    static void readFITSImage(fitsfile *fptr, NumMatrix<T>& M) {
    int naxis, status = 0;
    fits_get_img_dim(fptr, &naxis, &status);
    if (naxis!=2)
      throw std::invalid_argument("IO: naxis != 2. This is not a FITS image!");

    long naxes[2] = {1,1};
    fits_get_img_size(fptr, naxis, naxes, &status);
    M.resize(naxes[1],naxes[0]);
    long firstpix[2] = {1,1};
    T val;
    int imageformat = getFITSImageFormat(val);
    int datatype = getFITSDataType(val);
    fits_read_pix(fptr, datatype, firstpix, naxes[0]*naxes[1], NULL, M.c_array(), NULL, &status);
    if (status != 0)
      throw std::runtime_error("IO: Cannot read FITS image!");
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
      throw std::invalid_argument("IO: naxis != 2. This is not a FITS image!");

    long naxes[2] = {1,1};
    fits_get_img_size(fptr, naxis, naxes, &status);
    im.grid.setSize(0,0,naxes[0],naxes[1]);
    im.resize(im.grid.size());
    long firstpix[2] = {1,1};
    T val;
    int datatype = getFITSDataType(val);
    fits_read_pix(fptr, datatype, firstpix, im.size(), NULL, im.c_array(), NULL, &status);
    if (status != 0)
      throw std::runtime_error("IO: Cannot read FITS image!");
  }

  template <class T>
    static void readFITSImage(fitsfile *fptr, Image<complex<T> >& im) {
    int naxis, status = 0;
    fits_get_img_dim(fptr, &naxis, &status);
    if (naxis!=2)
      throw std::invalid_argument("IO: naxis != 2. This is not a FITS image!");

    long naxes[2] = {1,1};
    fits_get_img_size(fptr, naxis, naxes, &status);
    im.grid.setSize(0,0,naxes[0],naxes[1]);
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
      im(i) += complex<T>(0,component(i));
    if (status != 0)
      throw std::runtime_error("IO: Cannot read FITS image!");
  }

  /// Read in keyword from FITS header.
  /// For std::string keywords, use readFITSKeywordString() instead.
  template <class T>
    static void readFITSKeyword(fitsfile *fptr, std::string key, T& val) {
    int status = 0;
    char* comment = NULL;
    fits_read_key (fptr,getFITSDataType(val), const_cast<char *>(key.c_str()),&val,comment, &status);
    if (status != 0)
      throw std::invalid_argument("IO: Cannot read FITS keyword " + key + "!");  }

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
  static void writePPMImage(std::string filename,std::string colorscheme, std::string scaling, data_t min, data_t max, const Grid& grid, const NumVector<data_t>& data);
  /// Create RGB representation of data.
  /// This function is usefull for manipulations of the data in RGB space. 
  /// See writePPMImage() for details.
  static void makeRGBImage(NumMatrix<unsigned int>& rgbImage, std::string colorscheme, std::string scaling, data_t min, data_t max, const Grid& grid, const NumVector<data_t>& data);
  /// Write RGBImage from makeRGBImage() to PPM file.
  static void writeRGB2PPMImage (std::string filename, const Grid& grid, const NumMatrix<unsigned int>& rgbImage);
  /// Add uniform noise from noisemean to noisemean+noiselimit.
  static void addUniformNoise(NumVector<data_t>& data, data_t noisemean, data_t noiselimit);
  /// Add gaussian noise.
  static void addGaussianNoise(NumVector<data_t>& data, data_t noisemean, data_t noisesigma);
  /// Add poissonian noise.
  /// The noise of the \f$i\f$th pixel is drawn from a gaussian distribution with
  /// \f$ \sigma_i = \sqrt{\mu_n + data_i}\f$.
  static void addPoissonianNoise(NumVector<data_t>& data, data_t noisemean);
  /// Convolve input with a 3x3 Gaussian.
  static void convolveGaussian(const NumVector<data_t>& input, NumVector<data_t>& result, int width,int height);
 private:
  static int makeColorMatrix(NumMatrix<unsigned int>& m, std::string colorscheme);
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
    int IO::getFITSDataType<complex<float> >(const complex<float>& entry) {
    return TCOMPLEX;
  }
  template<> inline
    int IO::getFITSDataType<complex<double> >(const complex<double>& entry) {
    return TDBLCOMPLEX;
  }
  template<> inline
    int IO::getFITSDataType<std::string>(const std::string& entry) {
    return TSTRING;
  }
  
} // end namespace
#endif
