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

/// set imageformat and datatype according to the cfitsio definitions
template <class T> 
void setFITSTypes(T entry, int& imageformat, int& datatype) {
  imageformat = DOUBLE_IMG;
  datatype = TDOUBLE;
  float tf;
  int ti;
  bool tb;
  unsigned char tc;
  if (typeid(entry) == typeid(tf)) {
    imageformat = FLOAT_IMG;
    datatype = TFLOAT;
  }
  if (typeid(entry) == typeid(ti)) {
    imageformat = SHORT_IMG;
    datatype = TINT;
  }
  if (typeid(entry) == typeid(tb) || typeid(entry) == typeid(tc)) {
    imageformat = BYTE_IMG;
    datatype = TBYTE;
  }
}

/// creates a single HDU FITS file from the data defined on given grid
/// \todo add support for fits header keywords of variable type
template <class T>
void writeFITSFile(std::string filename, const Grid& grid, const NumVector<T>& data, std::map<std::string, std::string> keywords = std::map<std::string, std::string>()) {
  int dim0 = (int) ceil ((grid.getStopPosition(0) - grid.getStartPosition(0))/grid.getStepsize(0))+1;
  int dim1 = (int) ceil ((grid.getStopPosition(1) - grid.getStartPosition(1))/grid.getStepsize(1))+1;
  long naxis = 2;      
  long naxes[2] = { dim0, dim1 };
  long npixels = dim0*dim1;

  int status = 0;
  fitsfile *outfptr;
  filename = "!"+filename; // overwrite existing file if necessary
  // open fits file
  fits_create_file(&outfptr,filename.c_str(), &status);
  
  // define image format and dataformat according to cfitsio definitions
  int imageformat, datatype;
  setFITSTypes(data(0),imageformat,datatype);
  // create pHDU
  fits_create_img(outfptr, imageformat, naxis, naxes, &status);
  // write pixel data
  long firstpix[2] = {1,1};
  fits_write_pix(outfptr,datatype,firstpix,npixels,const_cast<T *>(data.c_array()), &status);
  // insert keywords
  keywords["CREATOR"]= "ShapeLens++";
  for( std::map<std::string,std::string>::iterator iter = keywords.begin(); 
       iter != keywords.end(); iter++ )
    fits_update_key (outfptr, TSTRING,  const_cast<char *>((*iter).first.c_str()), const_cast<char *>((*iter).second.c_str()), "", &status);
  fits_close_file(outfptr, &status);
}

template <class T>
void writeFITSFile(std::string filename, const NumMatrix<T>& M, std::map<std::string,std::string> keywords = std::map<std::string, std::string>()) {
  Grid grid(0,M.getRows()-1,1, 0, M.getColumns()-1, 1);
  writeFITSFile(filename,grid,M.vectorize(0),keywords);
}

/// adds extensions to an existing FITS file from the data on the given grid
template <class T>
void addFITSExtension(std::string filename, std::string extname, const Grid& grid, const NumVector<T>& data, std::map<std::string, std::string> keywords = std::map<std::string, std::string>()) {
    int dim0 = (int) ceil ((grid.getStopPosition(0) - grid.getStartPosition(0))/grid.getStepsize(0))+1;
  int dim1 = (int) ceil ((grid.getStopPosition(1) - grid.getStartPosition(1))/grid.getStepsize(1))+1;
  long naxis = 2;      
  long naxes[2] = { dim0, dim1 };
  long npixels = dim0*dim1;

  int status = 0;
  fitsfile *outfptr;
  fits_open_file(&outfptr, filename.c_str(), READWRITE, &status);
  // define image format and dataformat according to cfitsio definitions
  int imageformat, datatype;
  setFITSTypes(data(0),imageformat,datatype);
  // create extHDU
  fits_create_img(outfptr, imageformat, naxis, naxes, &status);
  // write pixel data 
  long firstpix[2] = {1,1};
  fits_write_pix(outfptr,datatype,firstpix,npixels,const_cast<T* >(data.c_array()), &status);
  // insert keywords
  keywords["EXTNAME"]= extname;
  keywords["CREATOR"]= "ShapeLens++";
  for( std::map<std::string,std::string>::iterator iter = keywords.begin(); 
       iter != keywords.end(); iter++ )
    fits_update_key (outfptr, TSTRING,  const_cast<char *>((*iter).first.c_str()), const_cast<char *>((*iter).second.c_str()), "", &status);
  fits_close_file(outfptr, &status);
}

template <class T>
void addFITSExtension(std::string filename, std::string extname, const NumMatrix<T>& M, std::map<std::string,std::string> keywords = std::map<std::string, std::string>()) {
  Grid grid(0,M.getRows()-1,1, 0, M.getColumns()-1, 1);
  addFITSExtension(filename,extname,grid,M.vectorize(0),keywords);
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
/// \f$ \sigma_i = \sigma_n + \sqrt{data_i}\f$.
void addPoissonianNoise(NumVector<data_t>& data, data_t noisemean, data_t noisesigma);
/// Convolve input with a 3x3 Gaussian.
void convolveGaussian(const NumVector<data_t>& input, NumVector<data_t>& result, int width,int height);
#endif
