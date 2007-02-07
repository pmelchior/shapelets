#include <Grid.h>

namespace ublas = boost::numeric::ublas;

Grid::Grid() {
}

Grid::Grid(double instart, double instop, double instepsize) {
  start = stop = stepsize = ublas::vector<double>(1);
  start(0) = instart;
  stop(0) = instop;
  stepsize(0) = instepsize;
  int lines =  (int) ceil((stop(0)-start(0))/stepsize(0))+1;
  grid = ublas::matrix<double>(lines,1);
  double x = start(0);
  for(int i =0; i<lines; i++) {
    grid(i,0) = x;
    x+=stepsize(0);
  }
}

Grid::Grid(double start0, double stop0, double stepsize0, double start1, double stop1, double stepsize1) {
  start = stop = stepsize = ublas::vector<double>(2);
  start(0) = start0;
  start(1) = start1;
  stop(0) = stop0;
  stop(1) = stop1;
  stepsize(0) = stepsize0;
  stepsize(1) = stepsize1;
  int number0 = (int) ceil((stop(0)-start(0))/stepsize(0))+1;
  int number1 = (int) ceil((stop(1)-start(1))/stepsize(1))+1;
  grid = ublas::matrix<double>(number0*number1,2);
  double x0 = start(0), x1 = start(1);
  int lines =0;
  // since FITS images are counted along x axis
  // the grid's first number is running, then the seconds
  for(int i=0;i<number1;i++) {
    for (int j =0; j<number0; j++) {
      grid(lines,0) = x0;
      grid(lines,1) = x1;
      //std::cout << x0 << " " << x1 << std::endl;
      x0+=stepsize(0);
      lines++;
    }
    x0 = start(0);
    x1+=stepsize(1);
  }
}

Grid& Grid::operator= (const Grid& h) {
  grid = h.grid;
  start = h.start;
  stop = h.stop;
  stepsize = h.stepsize;
  return *this;
}

double& Grid::operator() (const unsigned int i, const unsigned int j) {
  return grid(i,j);
}

Point2D& Grid::operator() (const unsigned int i) {
    ublas::matrix_row<ublas::matrix<double> > row (grid, i);
    point = row; 
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

int Grid::getSize(bool direction) const {
  return (int) ceil((stop(direction)-start(direction))/stepsize(direction))+1;
}

int Grid::size() const {
  return grid.size1();
}
