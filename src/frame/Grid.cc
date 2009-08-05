#include "../../include/frame/Grid.h"

using namespace shapelens;

// set up default (null) transformations
// for external calls to WCS via Grid
WCS::WCS() {
  CT = new NullTransformation<data_t>();
  CT_ = new NullTransformation<data_t>();
}

WCS::~WCS() {
  if (CT!=NULL)
    delete CT;
  if (CT_!=NULL)
    delete CT_;
}

WCS::WCS(const WCS& wcs) {
  CT = wcs.CT->clone();
  CT_ = wcs.CT_->clone();
}

WCS& WCS::operator=(const WCS& wcs) {
  if (CT!=NULL)
    delete CT;
  CT = wcs.CT->clone();
  if (CT_!=NULL)
    delete CT_;
  CT_ = wcs.CT_->clone();
}

void WCS::set(const CoordinateTransformation<data_t>& C) {
  if (CT!=NULL)
    delete CT;
  CT = C.clone();
  if (CT_!=NULL)
    delete CT_;
  CT_ = C.getInverse();
}

void WCS::reset() {
  if (CT!=NULL)
    delete CT;
  if (CT_!=NULL)
    delete CT_;
}

Grid::Grid() :
  N0(0),
  N1(0),
  wcs_set(false)
{
}

Grid::Grid(int start0, int start1, unsigned int N0, unsigned int N1) :
  start0(start0),
  start1(start1),
  N0(N0),
  N1(N1),
  wcs_set(false)
{
}

void Grid::setSize(int start0_, int start1_, unsigned int N0_, unsigned int N1_) {
  start0 = start0_;
  start1 = start1_;
  N0 = N0_;
  N1 = N1_; 
}

// this is the most often used access mode, so make it fast
data_t Grid::operator() (unsigned long i, bool direction) const {
  if (wcs_set) {
    Point<data_t> P(start0 + i%N0, start1 + i/N0);
    wcs.CT->transform(P);
    return P(direction);
  } else {
    long offset;
    if (direction) {
      offset = i/N0;
      return data_t(start1 + offset);
    }
    else {
      offset = i%N0;
      return data_t(start0 + offset);
    }
  }
}

Point<data_t> Grid::operator() (unsigned long i) const {
  if (wcs_set) {
    Point<data_t> P(start0 + i%N0, start1 + i/N0);
    wcs.CT->transform(P);
    return P;
  }
  else
    return Point<data_t>(operator()(i,0),operator()(i,1));
}

void Grid::setWCS(const CoordinateTransformation<data_t>& C) {
  wcs.set(C);
  wcs_set = true;
}

void Grid::resetWCS() {
  wcs.reset();
  wcs_set = false;
}

const WCS& Grid::getWCS() const {
  return wcs;
}

Rectangle<data_t> Grid::getBoundingBox() const {
  Rectangle<data_t> bbox;
  if (wcs_set) {
    bbox.ll = bbox.tr = operator()(0);
    Point<int> IC;
    Point<data_t> P;
    // for regular transformations, we only have to consider edgepoints
    for (IC(0) = start0; IC(0) < start0 + N0; IC(0)++) {
      IC(1) = start1;
      P = operator()(getPixel(IC)); // is always within image
      if (P(0)<bbox.ll(0))
	bbox.ll(0) = P(0);
      if (P(1)<bbox.ll(1))
	bbox.ll(1) = P(1);
      if (P(0)>bbox.tr(0))
	bbox.tr(0) = P(0);
      if (P(1)>bbox.tr(1))
	bbox.tr(1) = P(1);
      IC(1) = start1 + N1 - 1;
      P = operator()(getPixel(IC)); // is always within image
      if (P(0)<bbox.ll(0))
	bbox.ll(0) = P(0);
      if (P(1)<bbox.ll(1))
	bbox.ll(1) = P(1);
      if (P(0)>bbox.tr(0))
	bbox.tr(0) = P(0);
      if (P(1)>bbox.tr(1))
	bbox.tr(1) = P(1);
    }
    for (IC(1) = start1; IC(1) < start1 + N1; IC(1)++) {
      IC(0) = start0;
      P = operator()(getPixel(IC)); // is always within image
      if (P(0)<bbox.ll(0))
	bbox.ll(0) = P(0);
      if (P(1)<bbox.ll(1))
	bbox.ll(1) = P(1);
      if (P(0)>bbox.tr(0))
	bbox.tr(0) = P(0);
      if (P(1)>bbox.tr(1))
	bbox.tr(1) = P(1);
      IC(0) = start0 + N0 - 1;
      P = operator()(getPixel(IC)); // is always within image
      if (P(0)<bbox.ll(0))
	bbox.ll(0) = P(0);
      if (P(1)<bbox.ll(1))
	bbox.ll(1) = P(1);
      if (P(0)>bbox.tr(0))
	bbox.tr(0) = P(0);
      if (P(1)>bbox.tr(1))
	bbox.tr(1) = P(1);
    }
  }
  else {  
    bbox.ll = operator()(0);
    bbox.tr = operator()(Grid::size()-1);
  }
  return bbox;
}

int Grid::getStartPosition(bool direction) const {
  if (direction)
    return start1;
  else
    return start0;
}

int Grid::getStopPosition(bool direction) const {
  if (direction)
    return start1 + N1;
  else
    return start0 + N0;
}

unsigned int Grid::getSize(bool direction) const {
  if (direction)
    return N1;
  else
    return N0;
}

unsigned long Grid::size() const {
  return N0*N1;
}

Point<int> Grid::getCoords(unsigned long pixel) const {
  return Point<int>(start0 + pixel%N0,start1 + pixel/N0);
}

Point<int> Grid::getCoords(const Point<data_t>& P) const {
  if (wcs_set) {
    Point<data_t> P_ = P;
    wcs.CT_->transform(P_);
    return Point<int>((int)floor(P_(0)),(int)floor(P_(1)));
  } else
    return Point<int>((int)floor(P(0)),(int)floor(P(1)));
}

long Grid::getPixel(const Point<int>& P) const {
  if (P(0) >= start0 && P(0) < start0 + N0 && P(1) >= start1 && P(1) < start1 + N1)
    return (long) (P(0)-start0) + (P(1)-start1)*N0;
  else
    return -1;
}

long Grid::getNeighborPixel(const Point<int>& P, unsigned char direction) const {
  long index;
  int x = P(0), y = P(1);
  switch(direction) {
  case 0: 
    // the pixel itself
    index = y*N0 + x;
    break;
  case 1: 
    if (y<N1-1) index = (y+1)*N0 + x ;  // top
    else index = -1;
    break;
  case 2:
    if (y<N1-1 && x<N0-1) index = (y+1)*N0 + x + 1;  // top right
    else index = -1;
    break;
  case 3:
    if (x<N0-1) index = y*N0 + x + 1;  // right neighbour
    else index = -1;
    break;
  case 4: 
    if (y>0 && x<N0-1) index = (y-1)*N0 + x + 1;  // bottom right
    else index = -1;
    break;  
  case 5: 
    if (y>0) index = (y-1)*N0 + x;  // bottom
    else index = -1;
    break;
  case 6: 
    if (y>0 && x>0) index = (y-1)*N0 + x - 1;  // bottom left
    else index = -1;
    break;   
  case 7: 
    if (x>0) index = y*N0 + x - 1; // left
    else index = -1;
    break;
  case 8: 
    if (y<N1-1 && x>0) index = (y+1)*N0 + x - 1;  // top left
    else index = -1;
    break;  
  }
  return index;
}

long Grid::getNeighborPixel(unsigned long pixel, unsigned char direction) const {
  if (direction==0)
    return pixel;
  else
    return getNeighborPixel(getCoords(pixel),direction);
}

