#include <frame/CorrelationFunction.h>
#include <gsl/gsl_math.h>

typedef unsigned int uint;

CorrelationFunction::CorrelationFunction () {
}

CorrelationFunction::CorrelationFunction (const Image<data_t>& im, const SegmentationMap& segMap, uint insize, bool mask) {
  size = insize;
  uint x,y,x1,y1,i;

  std::map<uint, uint> index;
  makeIndexSetDistances(index);
  uint length = index.size();

  xi.resize(length);
  sigma.resize(length); // dist has been resized by makeIndexSetDistances()
  NumVector<uint> num(length);


  int axsize0 = im.getSize(0), axsize1 = im.getSize(1);
  for (i =0; i < im.size(); i++) {
    // choose only noise pixels
    if (!mask || segMap(i) <= 0) {
      im.getGrid().getCoords(i,x,y);
      for (x1=x-size; x1<= x+size; x1++) {
	for (y1=y-size; y1<= y+size; y1++) {
	  if (x1>=0 && x1<axsize0 && y1>=0 && y1 < axsize1) {
	    uint j = im.getGrid().getPixel(x1,y1);
	    // again: choose only noise pixels
	    if (!mask || segMap(j) <= 0) {
	      uint dist = (x1-x)*(x1-x)+(y1-y)*(y1-y);
	      xi(index[dist]) += im(i)*im(j);
	      num(index[dist])++;
	    }
	  }
	}
      }
    }
  }
  // normalize values by numbers of pixels counted per entry
  for (i=0; i<length; i++) {
    if (num(i) > 0)
      xi(i) /= num(i);
  }

  // compute the variance
  for (i =0; i < im.size(); i++) {
    // choose only noise pixels
    if (!mask || segMap(i) <= 0) {
      im.getGrid().getCoords(i,x,y);
      for (x1=x-size; x1<= x+size; x1++) {
	for (y1=y-size; y1<= y+size; y1++) {
	  if (x1>=0 && x1<axsize0 && y1>=0 && y1 < axsize1) {
	    uint j = im.getGrid().getPixel(x1,y1);
	    // again: choose only noise pixels
	    if (!mask || segMap(j) <= 0) {
	      uint dist = (x1-x)*(x1-x)+(y1-y)*(y1-y);
	      sigma(index[dist]) += gsl_pow_2(xi(index[dist]) - im(i)*im(j));
	    }
	  }
	}
      }
    }
  }
  // normalize values by numbers of pixels counted per entry
  for (i=0; i<length; i++) {
    if (num(i) > 0)
      sigma(i) = sqrt(sigma(i)/(num(i)*(num(i) - 1)));
  }
}

// same as above but with out segmentation map;
// this runs over all data pixels.
CorrelationFunction::CorrelationFunction (const NumVector<data_t>& data, const Grid& grid, uint insize) {
  size = insize;
  uint x,y,x1,y1,i;

  std::map<uint, uint> index;
  makeIndexSetDistances(index);
  uint length = index.size();

  xi.resize(length);
  sigma.resize(length); // dist has been resized by makeIndexSetDistances()
  NumVector<uint> num(length);

  int axsize0 = grid.getSize(0), axsize1 = grid.getSize(1);
  for (i =0; i < grid.size(); i++) {
    grid.getCoords(i,x,y);
    for (x1=x-size; x1<= x+size; x1++) {
      for (y1=y-size; y1<= y+size; y1++) {
	if (x1>=0 && x1<axsize0 && y1>=0 && y1 < axsize1) {
	  uint j = grid.getPixel(x1,y1);
	  uint dist = (x1-x)*(x1-x)+(y1-y)*(y1-y);
	  xi(index[dist]) += data(i)*data(j);
	  num(index[dist])++;
	}
      }
    }
  }
  // normalize values by numbers of pixels counted per entry
  for (i=0; i<length; i++) {
    if (num(i) > 0)
      xi(i) /= num(i);
  }

  // compute the variance
  for (i =0; i < grid.size(); i++) {
    grid.getCoords(i,x,y);
    for (x1=x-size; x1<= x+size; x1++) {
      for (y1=y-size; y1<= y+size; y1++) {
	if (x1>=0 && x1<axsize0 && y1>=0 && y1 < axsize1) {
	  uint j = grid.getPixel(x1,y1);
	  uint dist = (x1-x)*(x1-x)+(y1-y)*(y1-y);
	  sigma(index[dist]) += gsl_pow_2(xi(index[dist]) - data(i)*data(j));
	}
      }
    }
  }
  // normalize values by numbers of pixels counted per entry
  for (i=0; i<length; i++) {
    if (num(i) > 0)
      sigma(i) = sqrt(sigma(i)/(num(i)*(num(i) - 1)));
  }
}

void CorrelationFunction::operator= (const CorrelationFunction& xi2) {
  xi = xi2.xi;
  dist = xi2.dist;
  sigma = xi2.sigma;
  size = xi2.size;
}

const NumVector<data_t>& CorrelationFunction::getCorrelationFunction() const {
  return xi;
}
NumVector<data_t>& CorrelationFunction::accessCorrelationFunction() {
  return xi;
}

const NumVector<data_t>& CorrelationFunction::getDistances() const {
  return dist;
}

NumVector<data_t>& CorrelationFunction::accessDistances() {
  return dist;
}

const NumVector<data_t>& CorrelationFunction::getCorrelationError() const {
  return sigma;
}

NumVector<data_t>& CorrelationFunction::accessCorrelationError() {
  return sigma;
}

// setup map distance^2 -> vector index
// both distances and indices are ordered
void CorrelationFunction::makeIndexSetDistances(std::map<unsigned int, unsigned int>& index) {
  // first setup distances: creates a ordered list of keys, possible
  // distances between two pixels (smaller than rmax)
  for (uint i=0; i<=size; i++)
    for (uint j=i; j<=size; j++)
      index[i*i + j*j] = 0;
  // set keys to sequences of vector indices between 0 and length
  uint n=0;
  dist.resize(index.size());
  for( std::map<uint,uint>::iterator iter = index.begin(); iter != index.end(); iter++ ) {
    (*iter).second = n;
    dist(n) = sqrt((*iter).first);
    n++;
  }   
}
  
