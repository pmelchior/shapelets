#include <FitsImage.h>

using namespace std;

FitsImage::FitsImage(string infilename) {
  filename = infilename;
  read();
}

void FitsImage::setFilename(string infilename){
  filename = infilename;
  read();
}

std::string FitsImage::getFilename() {
  return filename;
}

void FitsImage::read() {
  fitsfile *fptr;
  int status = 0;
  fits_open_file(&fptr, filename.c_str(), READONLY, &status);
  int naxis;
  fits_get_img_dim(fptr, &naxis, &status);
  if (naxis!=2) {
    cout << "FitsImage: naxis != 2. This is not a FITS image!" << endl;
    terminate();
  } else {
    long naxes[2] = {1,1};
    fits_get_img_size(fptr, naxis, naxes, &status);
    axsize0 = naxes[0];
    axsize1 = naxes[1];
    grid = Grid(0,axsize0-1,1,0,axsize1-1,1);
    long npixels = axsize0*axsize1;
    data.resize(npixels);
    long firstpix[2] = {1,1};
    fits_read_pix(fptr, TDOUBLE, firstpix, npixels, NULL,data.c_array(), NULL, &status);
    fits_close_file(fptr, &status);
  }
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

void FitsImage::getCoords(uint pixel, int& x, int& y) {
  x = pixel%axsize0;
  y = pixel/axsize0;
}

unsigned int FitsImage::getPixel(int x, int y) {
  return (unsigned int) x + y*axsize0;
}
