#ifndef IO_H
#define IO_H

#include <fitsio.h>
#include <iostream>
#include <string>
#include <map>
#include <NumMatrix.h>
#include <Typedef.h>
#include <frame/Grid.h>

namespace shapelens {
template <class T>
class Image;

/// Functions for reading and writing into several formats
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
  template <class T> 
    static  int getFITSImageFormat(T entry) {
    char tc;
    unsigned char tuc;
    int ti;
    unsigned int tui;
    long tl;
    unsigned long tul;
    float tf;
    double td;
    if (typeid(entry) == typeid(tc) || typeid(entry) == typeid(tuc))
      return BYTE_IMG;
    if (typeid(entry) == typeid(ti))
      return SHORT_IMG;
    if (typeid(entry) == typeid(tui))
      return USHORT_IMG;
    if (typeid(entry) == typeid(tl))
      return LONG_IMG;
    if (typeid(entry) == typeid(tul))
      return ULONG_IMG;
    if (typeid(entry) == typeid(tf))
      return FLOAT_IMG;
    if (typeid(entry) == typeid(td))
      return DOUBLE_IMG;
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
  template <class T>
    static int getFITSDataType(T entry) {
    bool tb;
    char tc;
    unsigned char tuc;
    int ti;
    unsigned int tui;
    long tl;
    unsigned long tul;
    float tf;
    double td;
    complex<float> tcf;
    complex<double> tcd;
    std::string ts;
    if (typeid(entry) == typeid(tb) || typeid(entry) == typeid(tc) || typeid(entry) == typeid(tuc))
      return TBYTE;
    if (typeid(entry) == typeid(ti))
      return TINT;
    if (typeid(entry) == typeid(tui))
      return TUINT;
    if (typeid(entry) == typeid(tl))
      return TLONG;
    if (typeid(entry) == typeid(tul))
      return TULONG;
    if (typeid(entry) == typeid(tf))
      return TFLOAT;
    if (typeid(entry) == typeid(td))
      return TDOUBLE;
    if (typeid(entry) == typeid(tcf))
      return TCOMPLEX;
    if (typeid(entry) == typeid(tcd))
      return TDBLCOMPLEX;
    if (typeid(entry) == typeid(ts))
      return TSTRING;
  }

  /// Open FITS file.
  /// If <tt>write == false</tt>, the file will be opened in read-only mode.
  static fitsfile* openFITSFile(std::string filename, bool write=false);
  /// Create new FITS file.
  /// If the file \p filename already exists, it will be overwritten.
  static fitsfile* createFITSFile(std::string filename);
  /// Close FITS file pointer.
  static int closeFITSFile(fitsfile* fptr);
  /// Set/update std::string keyword in FITS file header.
  static int updateFITSKeywordString(fitsfile *outfptr, std::string keyword, std::string value, std::string comment="");
  /// Append \p history to FITS header histroy.
  static int appendFITSHistory(fitsfile *outfptr, std::string history);
  /// Read std::string keyword from FITS header.
  static int readFITSKeywordString(fitsfile *fptr, std::string key, std::string& val);
  /// Read FITS keyword cards directly.
  static int readFITSKeyCards(fitsfile *fptr, std::string key, std::string& value);

  /// Write FITS image from a NumVector and an appropriate Grid.
  /// The datatype will be automatically adjusted, based on the
  /// result of getFITSImageFormat() and getFITSDataType().
  /// If \p extname is non-empty, the keyword \p EXTNAME is set to \p extname.
  template <class T>
    static int writeFITSImage(fitsfile *outfptr, const Grid& grid, const NumVector<T>& data, std::string extname="") {
    int dim0 = grid.getSize(0);
    int dim1 = grid.getSize(1);
    long naxis = 2;      
    long naxes[2] = { dim0, dim1 };
    long npixels = dim0*dim1;

    // define image format and dataformat according to cfitsio definitions
    int imageformat = getFITSImageFormat(data(0));
    int datatype = getFITSDataType(data(0));
    // create HDU
    int status = 0;
    fits_create_img(outfptr, imageformat, naxis, naxes, &status);
    // write pixel data
    long firstpix[2] = {1,1};
    fits_write_pix(outfptr,datatype,firstpix,npixels,const_cast<T *>(data.c_array()), &status);
    // insert creator and extname keywords
    if (extname != "")
      status = updateFITSKeywordString (outfptr, "EXTNAME", extname);
    status = updateFITSKeywordString (outfptr, "CREATOR", "ShapeLens++");
    return status;
  }

  /// Write FITS image from an Image<T>.
  /// The datatype will be automatically adjusted, based on the
  /// result of getFITSImageFormat() and getFITSDataType().
  /// If \p extname is non-empty, the keyword \p EXTNAME is set to \p extname.
  template <class T>
    static int writeFITSImage(fitsfile *outfptr, const Image<T>& image, std::string extname="") {
    return writeFITSImage(outfptr,image.grid,image,extname);
  }

  /// Write FITS image from an NumMatrix<T>.
  /// The datatype will be automatically adjusted, based on the
  /// result of getFITSImageFormat() and getFITSDataType().
  /// If \p extname is non-empty, the keyword \p EXTNAME is set to \p extname.
  template <class T>
    static int writeFITSImage(fitsfile *outfptr, const NumMatrix<T>& M, std::string extname="") {
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
      status = updateFITSKeywordString (outfptr, "EXTNAME", extname);
    status = updateFITSKeywordString (outfptr, "CREATOR", "ShapeLens++");
    return status;
  }

  /// Set/update keyword in FITS file header.
  /// For setting string keywords, use updateFITSKeywordString() instead.
  template <class T>
    static int updateFITSKeyword(fitsfile *outfptr, std::string keyword, T value, std::string comment="") {
    int status = 0;
    fits_write_key (outfptr, getFITSDataType(value), const_cast<char *>(keyword.c_str()), &value, const_cast<char *>(comment.c_str()), &status);
    return status;
  }

  /// Read FITS image into NumMatrix<T>.
  /// \p M is adjusted to hold the contents of the image; the image value are 
  /// automatically casted to the type \p T of \p M.
  template <class T>
    static int readFITSImage(fitsfile *fptr, NumMatrix<T>& M) {
    int naxis, status = 0;
    fits_get_img_dim(fptr, &naxis, &status);
    if (naxis!=2) {
      std::cerr << "IO: naxis != 2. This is not a FITS image!" << std::endl;
      std::terminate();
    }
    long naxes[2] = {1,1};
    fits_get_img_size(fptr, naxis, naxes, &status);
    M.resize(naxes[1],naxes[0]);
    long firstpix[2] = {1,1};
    T val;
    int imageformat = getFITSImageFormat(val);
    int datatype = getFITSDataType(val);
    fits_read_pix(fptr, datatype, firstpix, naxes[0]*naxes[1], NULL, M.c_array(), NULL, &status);
    return status;
  }

  /// Read FITS image into NumVector<T> and a Grid.
  /// \p v is adjusted to hold the contents of the image; the image value are 
  /// automatically casted to the type \p T of \p v.\n
  /// \p grid is set to Grid(0,0,N,M), where \p N and \p M are the row and
  /// column numbers of the FITS image.
  template <class T>
    static int readFITSImage(fitsfile *fptr, Grid& grid, NumVector<T>& v) {
    int naxis, status = 0;
    fits_get_img_dim(fptr, &naxis, &status);
    if (naxis!=2) {
      std::cerr << "IO: naxis != 2. This is not a FITS image!" << std::endl;
      std::terminate();
    }
    long naxes[2] = {1,1};
    fits_get_img_size(fptr, naxis, naxes, &status);
    grid = Grid(0,0,naxes[0],naxes[1]);
    v.resize(grid.size());
    long firstpix[2] = {1,1};
    T val;
    int datatype = getFITSDataType(val);
    fits_read_pix(fptr, datatype, firstpix, grid.size(), NULL, v.c_array(), NULL, &status);
    return status;
  }
  /// Read FITS image into Image<T>.
  /// \p im is adjusted to hold the contents of the image; the image value are 
  /// automatically casted to the type \p T of \p im.\n
  /// \p Image<T>::grid is set to Grid(0,0,N,M), where \p N and \p M 
  /// are the row and column numbers of the FITS image.
  template <class T>
    static int readFITSImage(fitsfile *fptr, Image<T>& im) {
    int naxis, status = 0;
    fits_get_img_dim(fptr, &naxis, &status);
    if (naxis!=2) {
      std::cerr << "IO: naxis != 2. This is not a FITS image!" << std::endl;
      std::terminate();
    }
    long naxes[2] = {1,1};
    fits_get_img_size(fptr, naxis, naxes, &status);
    im.grid = Grid(0,0,naxes[0],naxes[1]);
    im.resize(im.grid.size());
    long firstpix[2] = {1,1};
    T val;
    int datatype = getFITSDataType(val);
    fits_read_pix(fptr, datatype, firstpix, im.grid.size(), NULL, im.c_array(), NULL, &status);
    return status;
  }

  /// Read in keyword from FITS header.
  /// For std::string keywords, use readFITSKeywordString() instead.
  template <class T>
    static int readFITSKeyword(fitsfile *fptr, std::string key, T& val) {
    int status = 0;
    char* comment = NULL;
    fits_read_key (fptr,getFITSDataType(val), const_cast<char *>(key.c_str()),&val,comment, &status);
    return status;
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
} // end namespace
#endif
