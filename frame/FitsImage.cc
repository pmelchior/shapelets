#include <FitsImage.h>
#include <cmath>
#include <valarray>

using namespace std;
using namespace CCfits;

FitsImage::FitsImage() {
  isRead = 0;
  axsize0 = axsize1 = 0;
}

FitsImage::FitsImage(string infilename) {
  filename = infilename;
  isRead = 0;
  axsize0 = axsize1 = 0;
}

void FitsImage::setFilename(string infilename){
  filename = infilename;
}
void FitsImage::read(int ext) {
  pInfile = auto_ptr<FITS>(new FITS(filename,Read,true));
  // read primary (ext=0) or extension
  if (ext==0) {
    PHDU& image = pInfile->pHDU();
    getContent(image);
  } else {
    ExtHDU& image = pInfile->extension(ext);
    getContent(image);
  }
  isRead = 1;
  grid = Grid(0,axsize0-1,1,0,axsize1-1,1);

}

template <class C>
void FitsImage::getContent (C& image) {
  image.readAllKeys();
  int bitpix;
  image.readKey("BITPIX",bitpix);
  // FIXME
  axsize0 = (unsigned int)image.axis(0);  // columns
  axsize1 = (unsigned int)image.axis(1);  // lines
  // no explicit typecast, just assign elements to a valarray<double>
  valarray<double> contents;
  image.read(contents);
  // copy to NumVector<double>
  data = contents;
}

const NumVector<double>& FitsImage::getData() {
  return data;
}

NumVector<double>& FitsImage::accessData() {
  return data;
}

const Grid& FitsImage::getGrid() {
  return grid;
}

unsigned int FitsImage::getNumberOfPixels() {
  return axsize0*axsize1;
}

unsigned int FitsImage::getSize(bool direction) {
  if (direction == 0) return axsize0;
  else return axsize1;
}

void FitsImage::printHeader() {
  if (isRead)
  std::cout << pInfile->pHDU() << std::endl;
}

void FitsImage::getCoords(uint pixel, int& x, int& y) {
  x = pixel%axsize0;
  y = pixel/axsize0;
}

unsigned int FitsImage::getPixel(int x, int y) {
  return (unsigned int) x + y*axsize0;
}
