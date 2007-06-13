#include <Grid.h>

Grid::Grid() {
}

Grid::Grid(double start0, double stop0, double stepsize0) {
  start = Point2D(start0,0.);
  stop = Point2D(stop0,0.);
  stepsize = Point2D(stepsize0,0.);
  axsize0 = computeSize(0);
  axsize1 = 1;
}

Grid::Grid(double start0, double stop0, double stepsize0, double start1, double stop1, double stepsize1) {
  start = Point2D(start0,start1);
  stop = Point2D(stop0,stop1);
  stepsize = Point2D(stepsize0,stepsize1);
  axsize0 = computeSize(0);
  axsize1 = computeSize(1);
}

double Grid::operator() (unsigned int index, bool direction) const {
  uint offset;
  if (direction)
    offset = index/axsize0;
  else
    offset = index%axsize0;
  return start(direction) + offset*stepsize(direction);
}

Point2D Grid::operator() (unsigned int i) const {
  Point2D point(operator()(i,0),operator()(i,1));
  return point;
}

  
double Grid::getStepsize(bool direction) const {
  return stepsize(direction);
}

double Grid::getStartPosition(bool direction) const {
  return start(direction);
}

double Grid::getStopPosition(bool direction) const {
  return stop(direction);
}

unsigned int Grid::getSize(bool direction) const {
  if (direction)
    return axsize1;
  else
    return axsize0;
}

unsigned int Grid::size() const {
  return axsize0*axsize1;
}

unsigned int Grid::computeSize(bool direction) const {
  return (int) ceil((stop(direction)-start(direction))/stepsize(direction))+1;
}

void Grid::getCoords(unsigned int pixel, unsigned int& x, unsigned int& y) const {
  x = pixel%getSize(0);
  y = pixel/getSize(0);
}

unsigned int Grid::getPixel(unsigned int x, unsigned int y) const {
  return (unsigned int) x + y*getSize(0);
}

int Grid::getNeighborPixel(unsigned int pixel, unsigned int x, unsigned int y, unsigned int direction) const {
  int index;
  unsigned int axsize0 = getSize(0), axsize1= getSize(1);
  switch(direction) {
  case 0: 
    // the pixel itself
    index = pixel;
    break;
  case 1: 
    if (y<axsize1-1) index = (y+1)*axsize0 + x ;  // top
    else index = -1;
    break;
  case 2:
    if (y<axsize1-1 && x<axsize0-1) index = (y+1)*axsize0 + x + 1;  // top right
    else index = -1;
    break;
  case 3:
    if (x<axsize0-1) index = y*axsize0 + x + 1;  // right neighbour
    else index = -1;
    break;
  case 4: 
    if (y>0 && x<axsize0-1) index = (y-1)*axsize0 + x + 1;  // bottom right
    else index = -1;
    break;  
  case 5: 
    if (y>0) index = (y-1)*axsize0 + x;  // bottom
    else index = -1;
    break;
  case 6: 
    if (y>0 && x>0) index = (y-1)*axsize0 + x - 1;  // bottom left
    else index = -1;
    break;   
  case 7: 
    if (x>0) index = y*axsize0 + x - 1; // left
    else index = -1;
    break;
  case 8: 
    if (y<axsize1-1 && x>0) index = (y+1)*axsize0 + x - 1;  // top left
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
