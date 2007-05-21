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
