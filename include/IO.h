#ifndef IO_H
#define IO_H

#include <fitsio.h>
#include <iostream>
#include <string>
#include <map>
#include <NumVector.h>
#include <NumMatrix.h>
#include <Typedef.h>
#include <frame/Grid.h>

/// Functions for reading and writing into several formats

template <class T> 
  int getFITSImageFormat(T entry) {
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
  


template <class T>
  int getFITSDataType(T entry) {
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

fitsfile* openFITSFile(std::string filename, bool write=0);
fitsfile* createFITSFile(std::string filename);
int closeFITSFile(fitsfile* fptr);
int updateFITSKeywordString(fitsfile *outfptr, std::string keyword, std::string value, std::string comment="");
int appendFITSHistory(fitsfile *outfptr, std::string history);
int readFITSKeywordString(fitsfile *fptr, std::string key, std::string& val);
int readFITSKeyCards(fitsfile *fptr, std::string key, std::string& value);

template <class T>
int writeFITSImage(fitsfile *outfptr, const Grid& grid, const NumVector<T>& data, std::string extname="") {
  int dim0 = (int) ceil ((grid.getStopPosition(0) - grid.getStartPosition(0))/grid.getStepsize(0))+1;
  int dim1 = (int) ceil ((grid.getStopPosition(1) - grid.getStartPosition(1))/grid.getStepsize(1))+1;
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

template <class T>
int writeFITSImage(fitsfile *outfptr, const NumMatrix<T>& M, std::string extname="") {
  Grid grid(0,M.getRows()-1,1, 0, M.getColumns()-1, 1);
  return writeFITSImage(outfptr,grid,M.vectorize(1),extname);
}

template <class T>
int updateFITSKeyword(fitsfile *outfptr, std::string keyword, T value, std::string comment="") {
  int status = 0;
  fits_write_key (outfptr, getFITSDataType(value), const_cast<char *>(keyword.c_str()), &value, const_cast<char *>(comment.c_str()), &status);
  return status;
}

template <class T>
int readFITSImage(fitsfile *fptr, NumMatrix<T>& M) {
  int naxis, status = 0;
  fits_get_img_dim(fptr, &naxis, &status);
  if (naxis!=2) {
    std::cerr << "IO: naxis != 2. This is not a FITS image!" << std::endl;
    std::terminate();
  }
  long naxes[2] = {1,1};
  fits_get_img_size(fptr, naxis, naxes, &status);
  M.resize(naxes[0],naxes[1]);
  long firstpix[2] = {1,1};
  T val;
  int imageformat = getFITSImageFormat(val);
  int datatype = getFITSDataType(val);
  fits_read_pix(fptr, datatype, firstpix, naxes[0]*naxes[1], NULL, M.c_array(), NULL, &status);
  return status;
}

template <class T>
int readFITSImage(fitsfile *fptr, Grid& grid, NumVector<T>& v) {
  int naxis, status = 0;
  fits_get_img_dim(fptr, &naxis, &status);
  if (naxis!=2) {
    std::cerr << "IO: naxis != 2. This is not a FITS image!" << std::endl;
    std::terminate();
  }
  long naxes[2] = {1,1};
  fits_get_img_size(fptr, naxis, naxes, &status);
  grid = Grid(0,naxes[0]-1,1,0,naxes[1]-1,1);
  v.resize(grid.size());
  long firstpix[2] = {1,1};
  T val;
  int datatype = getFITSDataType(val);
  fits_read_pix(fptr, datatype, firstpix, grid.size(), NULL, v.c_array(), NULL, &status);
  return status;
}

template <class T>
int readFITSKeyword(fitsfile *fptr, std::string key, T& val) {
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
/// min and max indicate the ends of the accepted range of values, values smaller
/// (larger) than min (max) are set to min (max).
void writePPMImage(std::string filename,std::string colorscheme, std::string scaling, data_t min, data_t max, const Grid& grid, const NumVector<data_t>& data);
/// Create RGB representation of data.
/// This function is usefull for manipulations of the data in RGB space. 
/// See writePPMImage() for details.
void makeRGBImage(NumMatrix<unsigned int>& rgbImage, std::string colorscheme, std::string scaling, data_t min, data_t max, const Grid& grid, const NumVector<data_t>& data);
/// Write RGBImage from makeRGBImage() to PPM file.
void writeRGB2PPMImage (std::string filename, const Grid& grid, const NumMatrix<unsigned int>& rgbImage);
/// Add uniform noise from noisemean to noisemean+noiselimit.
void addUniformNoise(NumVector<data_t>& data, data_t noisemean, data_t noiselimit);
/// Add gaussian noise.
void addGaussianNoise(NumVector<data_t>& data, data_t noisemean, data_t noisesigma);
/// Add poissonian noise.
/// The noise of the \f$i\f$th pixel is drawn from a gaussian distribution with
/// \f$ \sigma_i = \sqrt{\mu_n + data_i}\f$.
void addPoissonianNoise(NumVector<data_t>& data, data_t noisemean);
/// Convolve input with a 3x3 Gaussian.
void convolveGaussian(const NumVector<data_t>& input, NumVector<data_t>& result, int width,int height);

#endif
