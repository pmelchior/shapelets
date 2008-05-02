#include <frame/Grid.h>

Grid::Grid() :
  N0(0),
  N1(0)
{
}

Grid::Grid(data_t start0, data_t start1, int N0, int N1) : //, data_t stepsize0, data_t stepsize1) :
  start(start0,start1),
  stepsize(1,1),
  N0(N0),
  N1(N1),
  //stop(start0 + N0*stepsize0, start1 + N1*stepsize1)
  stop(start0 + N0 - 1, start1 + N1 - 1)
{
}

data_t Grid::operator() (unsigned int index, bool direction) const {
  int offset;
  if (direction)
    offset = index/N0;
  else
    offset = index%N0;
  return start(direction) + offset*stepsize(direction);
}

Point2D Grid::operator() (unsigned int i) const {
  return Point2D(operator()(i,0),operator()(i,1));
}

  
data_t Grid::getStepsize(bool direction) const {
  return stepsize(direction);
}

data_t Grid::getStartPosition(bool direction) const {
  return start(direction);
}

data_t Grid::getStopPosition(bool direction) const {
  return stop(direction);
}

unsigned int Grid::getSize(bool direction) const {
  if (direction)
    return N1;
  else
    return N0;
}

unsigned int Grid::size() const {
  return N0*N1;
}

void Grid::getCoords(unsigned int pixel, unsigned int& x, unsigned int& y) const {
  x = pixel%N0;
  y = pixel/N0;
}

unsigned int Grid::getPixel(unsigned int x, unsigned int y) const {
  return (unsigned int) x + y*N0;
}

int Grid::getNeighborPixel(unsigned int pixel, unsigned int x, unsigned int y, unsigned int direction) const {
  int index;
  switch(direction) {
  case 0: 
    // the pixel itself
    index = pixel;
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

int Grid::getNeighborPixel(unsigned int pixel, unsigned int direction) const {
  uint x,y;
  getCoords(pixel,x,y);
  return getNeighborPixel(pixel,x,y,direction);
}
