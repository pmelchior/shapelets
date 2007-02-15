#ifndef IO_H
#define IO_H

#include <fitsio.h>
#include <iostream>
#include <string>
#include <Grid.h>
#include <NumVector.h>
#include <NumMatrix.h>
#include <map>

/// Functions for reading and writing into several formats

/// set imageformat and datatype according to the cfitsio definitions
template <class T> 
void setFITSTypes(T entry, int& imageformat, int& datatype) {
  imageformat = DOUBLE_IMG;
  datatype = TDOUBLE;
  float tf;
  int ti;
  if (typeid(entry) == typeid(tf)) {
    imageformat = FLOAT_IMG;
    datatype = TFLOAT;
  }
  if (typeid(entry) == typeid(ti)) {
    imageformat = SHORT_IMG;
    datatype = TINT;
  }
}

/// creates a single HDU FITS file from the data defined on given grid
/// \todo add support for fits header keywords of variable type
template <class T>
void writeFITSFile(std::string filename, const Grid& grid, const NumVector<T>& data, std::map<std::string, std::string>& keywords) {
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
  keywords["CREATOR"]= "Shapelets++";
  for( std::map<std::string,std::string>::iterator iter = keywords.begin(); 
       iter != keywords.end(); iter++ )
    fits_update_key (outfptr, TSTRING,  const_cast<char *>((*iter).first.c_str()), const_cast<char *>((*iter).second.c_str()), "", &status);
  fits_close_file(outfptr, &status);
}

/// adds extensions to an existing FITS file from the data on the given grid
template <class T>
void addFITSExtension(std::string filename, std::string extname, const Grid& grid, const NumVector<T>& data, std::map<std::string, std::string>& keywords) {
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
  keywords["CREATOR"]= "Shapelets++";
  for( std::map<std::string,std::string>::iterator iter = keywords.begin(); 
       iter != keywords.end(); iter++ )
    fits_update_key (outfptr, TSTRING,  const_cast<char *>((*iter).first.c_str()), const_cast<char *>((*iter).second.c_str()), "", &status);
  fits_close_file(outfptr, &status);
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
void writePPMImage(std::string filename,std::string colorscheme, std::string scaling, double min, double max, const Grid& grid, const NumVector<double>& data);
/// Create RGB representation of data.
/// This function is usefull for manipulations of the data in RGB space. 
/// See writePPMImage() for details.
void makeRGBImage(NumMatrix<unsigned int>& rgbImage, std::string colorscheme, std::string scaling, double min, double max, const Grid& grid, const NumVector<double>& data);
/// Write RGBImage from makeRGBImage() to PPM file.
void writeRGB2PPMImage (std::string filename, const Grid& grid, const NumMatrix<unsigned int>& rgbImage);
/// Add uniform noise from noisemean to noisemean+noiselimit.
void addUniformNoise(NumVector<double>& data, double noisemean, double noiselimit);
/// Add gaussian noise.
void addGaussianNoise(NumVector<double>& data, double noisemean, double noisesigma);
/// Add poissonian noise.
/// The noise of the \f$i\f$th pixel is drawn from a gaussian distribution with
/// \f$ \sigma_i = \sigma_n + \sqrt{data_i}\f$.
void addPoissonianNoise(NumVector<double>& data, double noisemean, double noisesigma);
/// Convolve input with a 3x3 Gaussian.
void convolveGaussian(const NumVector<double>& input, NumVector<double>& result, int width,int height);
#endif
