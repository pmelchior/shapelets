#include <frame/CorrelationFunction.h>
#include <gsl/gsl_math.h>

typedef unsigned int uint;
using namespace shapelens;

CorrelationFunction::CorrelationFunction () {
}

CorrelationFunction::CorrelationFunction (const Image<data_t>& data, data_t threshold, int limit) {
  maxLength = 1;
  int sig_pixels;
  compute(data);
  sig_pixels = xi.size();
  applyThreshold(threshold); // cut entries in xi which are below threshold
  // continue until some entries go below threshold
  while (sig_pixels == xi.size()) {
    if ((limit > 0 && maxLength < limit) || limit <= 0) {
      maxLength++;
      compute(data);
      sig_pixels = xi.size();
      applyThreshold(threshold);
    } else 
      break;
  }
}

CorrelationFunction::CorrelationFunction (const Image<data_t>& data, const SegmentationMap& segMap, data_t threshold, int limit) {
  maxLength = 1;
  int sig_pixels;
  compute(data,segMap);
  sig_pixels = xi.size();
  applyThreshold(threshold); // cut entries in xi which are below threshold
  // continue until some entries go below threshold
  while (sig_pixels == xi.size()) {
    if ((limit > 0 && maxLength < limit) || limit <= 0) {
      maxLength++;
      compute(data,segMap);
      sig_pixels = xi.size();
      applyThreshold(threshold);
    } else 
      break;
  }
}


CorrelationFunction::CorrelationFunction(const NumMatrix<data_t>& corr)  {
  if (corr.getRows()%2 != 1 || corr.getRows() != corr.getColumns()) {
    std::cerr << "CorrelationFunction: correlation matrix must be square and have odd dimension" << std::endl;
    std::terminate();
  }
  maxLength = (corr.getColumns() - 1)/2;
  Point2D<grid_t> p;
  int i,j;
  for (int i=0; i < corr.getRows(); i++) {
    for (int j =0; j < corr.getColumns(); j++) {
      p(0) = i - maxLength;
      p(1) = j - maxLength;
      xi[p] = corr(i,j);
      sigma[p] = 0;
    }
  }
}

void CorrelationFunction::compute(const Image<data_t>& im, const SegmentationMap& segMap) {
  int x,y,x1,y1,i,j;
  Point2D<grid_t> p;
  setPoints();
  grid_t startx = im.grid.getStartPosition(0), starty = im.grid.getStartPosition(1), stopx = im.grid.getStopPosition(0), stopy = im.grid.getStopPosition(1);

  // 1) compute mean of correlation
  for (i =0; i < im.size(); i++) {
    // choose only noise pixels
    if (segMap(i) == 0) {
      im.grid.getCoords(i,x,y);
      for (x1=x-maxLength; x1<= x+maxLength; x1++) {
	for (y1=y-maxLength; y1<= y+maxLength; y1++) {
	  if (x1>=startx && x1<stopx && y1>=starty && y1 < stopy) {
	    j = im.grid.getPixel(x1,y1);
	    // again: choose only noise pixels
	    if (segMap(j) == 0) {
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
    if (segMap(i) == 0) {
      im.grid.getCoords(i,x,y);
      for (x1=x-maxLength; x1<= x+maxLength; x1++) {
	for (y1=y-maxLength; y1<= y+maxLength; y1++) {
	  if (x1>=startx && x1<stopx && y1>=starty && y1 < stopy) {
	    j = im.grid.getPixel(x1,y1);
	    // again: choose only noise pixels
	    if (segMap(j) == 0) {
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

void CorrelationFunction::compute (const Image<data_t>& im) {
  int x,y,x1,y1,i,j;
  Point2D<grid_t> p;
  setPoints();
  grid_t startx = im.grid.getStartPosition(0), starty = im.grid.getStartPosition(1), stopx = im.grid.getStopPosition(0), stopy = im.grid.getStopPosition(1);

  // 1) compute mean of correlation
  for (i =0; i < im.size(); i++) {
    im.grid.getCoords(i,x,y);
    for (x1=x-maxLength; x1<= x+maxLength; x1++) {
      for (y1=y-maxLength; y1<= y+maxLength; y1++) {
	if (x1 >= startx && x1 < stopx && y1 >= starty && y1 < stopy) {
	  j = im.grid.getPixel(x1,y1);
	  p(0) = x1-x;
	  p(1) = y1-y;
	  xi[p] += im(i)*im(j);
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
  for (i =0; i < im.size(); i++) {
    im.grid.getCoords(i,x,y);
    for (x1=x-maxLength; x1<= x+maxLength; x1++) {
      for (y1=y-maxLength; y1<= y+maxLength; y1++) {
	if (x1 >= startx && x1 < stopx && y1 >= starty && y1 < stopy) {
	  j = im.grid.getPixel(x1,y1);
	  p(0) = x1-x;
	  p(1) = y1-y;
	  sigma[p] += gsl_pow_2(xi[p] - im(i)*im(j));
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


NumMatrix<data_t> CorrelationFunction::getCorrelationMatrix() const {
  NumMatrix<data_t> corr(2*maxLength+1, 2*maxLength+1);
  std::map<Point2D<grid_t>, data_t>::const_iterator iter;
  int i,j;
  for (iter = xi.begin(); iter != xi.end(); iter++) {
    i = iter->first(0) + maxLength;
    j = iter->first(1) + maxLength;
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
  for (int i = -maxLength; i <= maxLength; i++) {
    for (int j = -maxLength; j <= maxLength; j++) {
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
  for (int i = -maxLength; i <= maxLength; i++) {
    for (int j = -maxLength; j <= maxLength; j++) {
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
  return maxLength;
}

unsigned int CorrelationFunction::getSize() const {
  return xi.size();
}
