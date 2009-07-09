#include "../../include/frame/Grid.h"

using namespace shapelens;

Grid::Grid() :
  N0(0),
  N1(0),
  pos_set(false)
{
}

Grid::Grid(int start0, int start1, unsigned int N0, unsigned int N1) :
  start0(start0),
  start1(start1),
  N0(N0),
  N1(N1),
  pos_set(false)
{
}

data_t Grid::operator() (unsigned long i, bool direction) const {
  if (pos_set) {
    return pos[i](direction);
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
  if (pos_set)
    return pos[i];
  else
    return Point<data_t>(operator()(i,0),operator()(i,1));
}

void Grid::apply(const CoordinateTransformation<data_t>& C) {
  Point<data_t> P;
  for (unsigned long i=0; i < Grid::size(); i++) {
    P = Grid::operator()(i);
    C.transform(P);
    pos.push_back(P);
  }
  pos_set = true;
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

unsigned long Grid::getPixel(const Point<int>& P) const {
  return (unsigned long) (P(0)-start0) + (P(1)-start1)*N0;
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

