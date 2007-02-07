#ifndef IO_H
#define IO_H

#include <CCfits/CCfits>
#include <iostream>
#include <string>
#include <Grid.h>
#include <NumVector.h>
#include <NumMatrix.h>
#include <map>

/// Function for reading and writing into several formats

/// creates a single HDU FITS file from the data defined on given grid
template <class T>
void writeFITSFile(std::string filename, const Grid& grid, const NumVector<T>& data, std::map<std::string, std::string>& keywords) {
  int dim0 = (int) ceil ((grid.getStopPosition(0) - grid.getStartPosition(0))/grid.getStepsize(0))+1;
  int dim1 = (int) ceil ((grid.getStopPosition(1) - grid.getStartPosition(1))/grid.getStepsize(1))+1;
  long naxis = 2;      
  long naxes[2] = { dim0, dim1 };
  std::auto_ptr<CCfits::FITS> pFits(0);
  // since the call to new FITS(filename...) doesn't overwrite the file
  // we have to remove it before
  std::ostringstream filename_remove;
  filename_remove << "!" << filename;
  // define image format
  int imageformat = DOUBLE_IMG;
  float tf;
  int ti;
  if (typeid(data(0)) == typeid(tf))
    imageformat = FLOAT_IMG;
  if (typeid(data(0)) == typeid(ti))
    imageformat = SHORT_IMG;
  pFits.reset(new CCfits::FITS(filename_remove.str(),imageformat,naxis,naxes));
  // insert kewords
  keywords["CREATOR"]= "Shapelets++";
  for( std::map<std::string,std::string>::iterator iter = keywords.begin(); 
       iter != keywords.end(); iter++ )
    pFits->pHDU().addKey((*iter).first,(*iter).second,"");

  // copy to valarray
  valarray<T> array(data.size());
  for (int i=0; i < data.size(); i++)
    array[i] = data(i);
  // writing into fits file
  pFits->pHDU().write(1,data.size(),array);
}

/// adds extensions to an existing FITS file from the data on the given grid
template <class T>
void addFITSExtension(std::string filename, std::string extname, const Grid& grid, const NumVector<T>& data, std::map<std::string, std::string>& keywords) {
    int dim0 = (int) ceil ((grid.getStopPosition(0) - grid.getStartPosition(0))/grid.getStepsize(0))+1;
  int dim1 = (int) ceil ((grid.getStopPosition(1) - grid.getStartPosition(1))/grid.getStepsize(1))+1;
  std::vector<long> naxes(2);
  naxes[0] = dim0;
  naxes[1] = dim1;
  // define image format
  int imageformat = DOUBLE_IMG;
  float tf;
  int ti;
  if (typeid(data(0)) == typeid(tf))
    imageformat = FLOAT_IMG;
  if (typeid(data(0)) == typeid(ti))
    imageformat = SHORT_IMG;
  // open the existing FITS file
  CCfits::FITS infile(filename, CCfits::Write);
  // select next extension
  CCfits::ExtHDU* imageExt = infile.addImage(extname,imageformat,naxes);
  // copy to valarray
  valarray<T> array(data.size());
  for (int i=0; i < data.size(); i++)
    array[i] = data(i);
  // writing into fits file
  imageExt->write(1,data.size(),array);
  // insert kewords
  keywords["CREATOR"]= "Shapelets++";
  for( std::map<std::string,std::string>::iterator iter = keywords.begin();
       iter != keywords.end(); iter++ )
    imageExt->addKey((*iter).first,(*iter).second,"");
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
