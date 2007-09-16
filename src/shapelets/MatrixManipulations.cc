#include <shapelets/MatrixManipulations.h>
#include <math.h>
#include <iostream>

namespace ublas = boost::numeric::ublas;
typedef complex<double> Complex;

// maps the range m = -n, -n +2, .. , n to
// m = 0,1,..,n to efficiently store in matrix
unsigned int mIndex(int m, int n) {
  return (m + n)/2;
}

