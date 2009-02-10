#include <frame/CorrelationFunction.h>
#include <gsl/gsl_math.h>

typedef unsigned int uint;

CorrelationFunction::CorrelationFunction () {
}

CorrelationFunction::CorrelationFunction (const Image<data_t>& im, const SegmentationMap& segMap, int size, bool mask) :
  size(size) {
  int x,y,x1,y1,i,j;
  Point2D<grid_t> p;
  setPoints();
  int axsize0 = im.getSize(0), axsize1 = im.getSize(1);

  // 1) compute mean of correlation
  for (i =0; i < im.size(); i++) {
    // choose only noise pixels
    if (!mask || segMap(i) == 0) {
      im.grid.getCoords(i,x,y);
      for (x1=x-size; x1<= x+size; x1++) {
	for (y1=y-size; y1<= y+size; y1++) {
	  if (x1>=0 && x1<axsize0 && y1>=0 && y1 < axsize1) {
	    j = im.grid.getPixel(x1,y1);
	    // again: choose only noise pixels
	    if (!mask || segMap(j) == 0) {
	      p(0) = x1-x;
	      p(1) = y1-y;
	      xi[p] += im(i)*im(j);
	      num[p]++;
	    }
	  }
	}
      }
    }
  }
  // normalize values by numbers of pixels counted per entry
  std::map<Point2D<grid_t>, data_t>::iterator iter;
  for (iter = xi.begin(); iter != xi.end(); iter++) {
    uint n = num[iter->first];
    if (n > 0)
      iter->second /= n;
  }

  // 2) compute std of correlation given mean from above
  for (i =0; i < im.size(); i++) {
    // choose only noise pixels
    if (!mask || segMap(i) == 0) {
      im.grid.getCoords(i,x,y);
      for (x1=x-size; x1<= x+size; x1++) {
	for (y1=y-size; y1<= y+size; y1++) {
	  if (x1>=0 && x1<axsize0 && y1>=0 && y1 < axsize1) {
	    j = im.grid.getPixel(x1,y1);
	    // again: choose only noise pixels
	    if (!mask || segMap(j) == 0) {
	      p(0) = x1-x;
	      p(1) = y1-y;
	      sigma[p] += gsl_pow_2(xi[p] - im(i)*im(j));
	    }
	  }
	}
      }
    }
  }
  // normalize values by numbers of pixels counted per entry
  for (iter = sigma.begin(); iter != sigma.end(); iter++) {
    uint n = num[iter->first];
    if (n > 0)
      iter->second = sqrt((iter->second)/(n*(n-1)));
  }
}

CorrelationFunction::CorrelationFunction (const NumVector<data_t>& data, const Grid& grid, int size) :
  size(size) {
  int x,y,x1,y1,i,j;
  Point2D<grid_t> p;
  setPoints();
  grid_t startx = grid.getStartPosition(0), starty = grid.getStartPosition(1),
    stopx = grid.getStopPosition(0), stopy = grid.getStopPosition(1);

  // 1) compute mean of correlation
  for (i =0; i < data.size(); i++) {
    grid.getCoords(i,x,y);
    for (x1=x-size; x1<= x+size; x1++) {
      for (y1=y-size; y1<= y+size; y1++) {
	if (x1 >= startx && x1 < stopx && y1 >= starty && y1 < stopy) {
	  j = grid.getPixel(x1,y1);
	  p(0) = x1-x;
	  p(1) = y1-y;
	  xi[p] += data(i)*data(j);
	  num[p]++;
	}
      }
    }
  }
  // normalize values by numbers of pixels counted per entry
  std::map<Point2D<grid_t>, data_t>::iterator iter;
  for (iter = xi.begin(); iter != xi.end(); iter++) {
    uint n = num[iter->first];
    //std::cout << iter->first << "\t" << iter->second << "\t" << n << std::endl;
    if (n > 0)
      iter->second /= n;
  }

  // 2) compute std of correlation given mean from above
  for (i =0; i < data.size(); i++) {
    grid.getCoords(i,x,y);
    for (x1=x-size; x1<= x+size; x1++) {
      for (y1=y-size; y1<= y+size; y1++) {
	if (x1 >= startx && x1 < stopx && y1 >= starty && y1 < stopy) {
	  j = grid.getPixel(x1,y1);
	  p(0) = x1-x;
	  p(1) = y1-y;
	  sigma[p] += gsl_pow_2(xi[p] - data(i)*data(j));
	}
      }
    }
  }
  // normalize values by numbers of pixels counted per entry
  for (iter = sigma.begin(); iter != sigma.end(); iter++) {
    uint n = num[iter->first];
    if (n > 0)
      iter->second = sqrt((iter->second)/(n*(n-1)));
  }
}

CorrelationFunction::CorrelationFunction(const NumMatrix<data_t>& corr)  {
  if (size%2 != 1) {
    std::cerr << "CorrelationFunction: correlation matrix must have odd dimension" << std::endl;
    std::terminate();
  }
  size = (corr.getColumns() - 1)/2;
  Point2D<grid_t> p;
  int i,j;
  for (int i=0; i < corr.getRows(); i++) {
    for (int j =0; j < corr.getColumns(); j++) {
      p(0) = i - size;
      p(1) = j - size;
      xi[p] = corr(i,j);
      sigma[p] = 0;
    }
  }
}

void CorrelationFunction::operator= (const CorrelationFunction& xi2) {
  xi = xi2.xi;
  sigma = xi2.sigma;
  size = xi2.size;
}

NumMatrix<data_t> CorrelationFunction::getCorrelationMatrix() const {
  NumMatrix<data_t> corr(2*size+1, 2*size+1);
  std::map<Point2D<grid_t>, data_t>::const_iterator iter;
  int i,j;
  for (iter = xi.begin(); iter != xi.end(); iter++) {
    i = iter->first(0) + size;
    j = iter->first(1) + size;
    corr(i,j) = iter->second;
  }
  return corr;
}

const std::map<Point2D<grid_t>, data_t>& CorrelationFunction::getCorrelationFunction() const {
  return xi;
}
const std::map<Point2D<grid_t>, data_t>& CorrelationFunction::getCorrelationError() const {
  return sigma;
}
  
void CorrelationFunction::setPoints() {
  Point2D<grid_t> p;
  for (int i = -size; i <= size; i++) {
    for (int j = -size; j <= size; j++) {
      p(0) = i;
      p(1) = j;
      xi[p] = sigma[p] = 0;
      num[p] = 0;
    }
  }
}

void CorrelationFunction::applyThreshold(data_t thresh) {
  Point2D<grid_t> p;
  // don't use iterators here because the erase messes them up
  for (int i = -size; i <= size; i++) {
    for (int j = -size; j <= size; j++) {
      p(0) = i;
      p(1) = j;
      if (xi[p] < thresh*sigma[p]) {
	xi.erase(p);
	sigma.erase(p);
      }
    }
  }
}

unsigned int CorrelationFunction::getMaxLength() const {
  return size;
}
