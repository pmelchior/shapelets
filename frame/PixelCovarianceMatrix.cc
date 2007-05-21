#include <PixelCovarianceMatrix.h>
#include <map>
#include <vector>
#include <gsl/gsl_math.h>

typedef unsigned int uint;

PixelCovarianceMatrix::PixelCovarianceMatrix() {
  bandwidth = 0;
}

PixelCovarianceMatrix::PixelCovarianceMatrix(unsigned int bandwidth) {
  setBandwidth(bandwidth);
}

uint PixelCovarianceMatrix::getBandwidth() const {
  return bandwidth;
}

void PixelCovarianceMatrix::setBandwidth(uint bw) {
  bandwidth = bw;
  offset = NumVector<int>(bandwidth);
  entry = NumVector<double>(bandwidth);
}

void PixelCovarianceMatrix::setBand(unsigned int band, int os, double val) {
  if (band<bandwidth) {
    offset(band) = os;
    entry(band) = val;
  }
}

double PixelCovarianceMatrix::operator()(unsigned int i, unsigned int j) const {
  bool found = 0;
  int k;
  for (uint b=0; b <bandwidth; b++) {
    k = i + offset(b);
    if (k==j) {
      found = 1;
      break;
    }
  }
  if (found)
    return entry(j);
  else
    return 0;
}

void PixelCovarianceMatrix::getCorrelationFunction(const NumVector<double>& data, const SegmentationMap& segMap, NumVector<double>& corr, NumVector<uint>& distance) const {
  uint x,y,x1,y1;
  uint size = 5; // the pixel 'radius' to consider
  uint length = 15; // the pyramid number of size

  NumVector<double> corr_mean(length);
  NumVector<int> corr_num(length);

  // setup map distance -> vector index
  std::map<uint, uint> index;
  uint n =0;
  for (uint i=0; i<size; i++) {
    for (uint j=i; j<size; j++) {
      index[i*i + j*j] = n;
      n++;
    }
  }

  //const SegmentationMap& segMap = obj.getSegmentationMap();
  //const NumVector<double>& data = obj.getData();
  int axsize0 = segMap.getSize(0), axsize1 = segMap.getSize(1);
  for (int i =0; i < data.size(); i++) {
    // choose only noise pixels
    if (segMap(i) <= 0) {
      segMap.getCoords(i,x,y);
      for (x1=x-size+1; x1< x+size; x1++) {
	for (y1=y-size+1; y1< y+size; y1++) {
	  if (x1>=0 && x1<axsize0 && y1>=0 && y1 < axsize1) {
	    uint j = segMap.getPixel(x1,y1);
	    // again: choose only noise pixels
	    if (segMap(j) <= 0) {
	      uint dist = (x1-x)*(x1-x)+(y1-y)*(y1-y);
	      corr_mean(index[dist]) += data(i)*data(j);
	      corr_num(index[dist])++;
	    }
	  }
	}
      }
    }
  }
  // normalize values by numbers of pixels counted per entry
  uint i;
  for (i=0; i<length; i++) {
    if (corr_num(i) > 0)
     corr_mean(i) /= corr_num(i);
  }
  
  // since map is ordered in distances not in vector indices
  // we have to turn things around
  corr.resize(length);
  distance.resize(length);
  i=0;
  for( std::map<uint,uint>::iterator iter = index.begin(); iter != index.end(); iter++ ) {
    distance(i) = (*iter).first;
    corr(i) = corr_mean((*iter).second);
    i++;
  }
}

void PixelCovarianceMatrix::setCovarianceMatrix(const NumVector<double>& data, const SegmentationMap& segMap, History& hist) {
  NumVector<double> corr;
  NumVector<uint> distance;
  hist.append("# Computing correlation function:\n");
  getCorrelationFunction(data,segMap,corr,distance);
  hist.append("# dist\tcorr\n");
  char diststr[100];
  char corrstr[100];
  for (uint i=0; i < corr.size(); i++) {
    sprintf(diststr,"%1.2f",sqrt(distance(i))); 
    sprintf(corrstr,"%.2f",corr(i));
    hist.append("# "+std::string(diststr)+"\t"+std::string(corrstr)+"\n");
  }

  // the bandwidth depends on the values of the correlation function
  // for distances larger than 0.
  // use 10% of the maximum (distance 0) value as threshold;
  // effectively zero pixel correlation
  if (corr(1) < 0.1*corr(0)) {
    bandwidth = 1;
    hist.append("# Restricting covariance matrix to 1 band: corr(r>0) = 0\n");
  }
  // larger correlation: 
  // consider now also top-right/top-left/bottom-right/bottom-left neighbor
  // correlation length of 2 pixels and larger should not happen
  // also for simplicity we restrict the cov. matrix to bandwidth of 9
  else {
    bandwidth = 9;
    hist.append("# Restricting covariance matrix to 9 bands: corr(r>=2) = 0\n");
  }
  setBandwidth(bandwidth);

  // set the bands of the covariance map according to bandwidth
  uint N1 = segMap.getSize(0);
  if (bandwidth ==1) {
    offset(0) = 0;
    entry(0) = corr(0);
  }
  if (bandwidth == 9) {
    offset(0) = -N1-1;
    offset(1) = -N1;
    offset(2) = -N1+1;
    offset(3) = -1;
    offset(4) = 0;
    offset(5) = 1;
    offset(6) = N1-1;
    offset(7) = N1;
    offset(8) = N1+1;
    entry(0) = entry(2) = entry(6) = entry(8) = corr(2);
    entry(1) = entry(3) = entry(5) = entry(7) = corr(1);
    entry(4) = corr(0);
  }
}

NumMatrix<double> PixelCovarianceMatrix::getMatrix(unsigned int N) const {
  NumMatrix<double> M(N,N);
  for (uint i=0; i<N; i++) {
    for (uint band=0; band<bandwidth; band++) {
      uint j = i + offset(band);
      if (j>=0 &&  j<N) {
	M(i,j) = entry(band);
      }
    }
  }
  return M;
}

PixelCovarianceMatrix PixelCovarianceMatrix::invert() const {
  // in principle we need to invert a (N1*N2)x(N1*N2) matrix here
  // since this entries on the bands are constant the size of the matrix does not matter
  // only entries and offsets are important
  // 1) set up small portion of true covariance matrix
  // 2) invert it using NumMatrix::invert();
  // 3) read off entries and offsets in the center of the matrix to avoid boundary effects
  //    set a threshold of 1% of the maximum value to be considered

  if (bandwidth==1) {
    PixelCovarianceMatrix inv(1);
    inv.setBand(0,0,1./entry(0));
    return inv;
  }
  else {
    uint N = GSL_MAX_INT(10,offset.max())*bandwidth;
    if (N%2==1) N++; // ensure that divides by 2
    NumMatrix<double> M = getMatrix(N).invert();
    writeFITSFile("V.fits",getMatrix(N));
    addFITSExtension("V.fits","INV",M);

    int center = N/2;
    double max = fabs(M(center,center));
    uint inv_bandwidth = 0;
    std::vector<int> inv_offset;
    std::vector<double> inv_entry;
    for (int j=1; j<N-1; j++) {
      double val = M(N-j,j);
      if (fabs(val) > 0.001*max) {
	inv_offset.push_back(2*(j-center));
	inv_entry.push_back(val);
	inv_bandwidth++;
      }
      val = M(N-j,j+1);
      if (fabs(val) > 0.001*max) {
	inv_offset.push_back(2*(j-center)+1);
	inv_entry.push_back(val);
	inv_bandwidth++;
      }
    }
    PixelCovarianceMatrix inv(inv_bandwidth);
    for (uint band=0; band<inv_bandwidth; band++)
      inv.setBand(band,inv_offset[band],inv_entry[band]);
    return inv;
  }
}

PixelCovarianceMatrix PixelCovarianceMatrix::transpose() const {
  PixelCovarianceMatrix trans(bandwidth);
  for (uint i=0; i<bandwidth; i++)
    trans.setBand(i,offset(bandwidth-i-1),entry(bandwidth-i-1));
  return trans;
}
  

NumVector<double> PixelCovarianceMatrix::operator*(const NumVector<double>& v) const {
  // diagonal matrix
  if (bandwidth==1) {
    NumVector<double> w=v;
    if (entry(0) != 1) 
      w *= entry(0);
    return w;
  } 
  // real banded matrix: do row-wise multiplication
  // exploit limited index range for inner loop
  else {
    int N = v.size();
    NumVector<double> w(N);
    for (int i=0; i<N; i++)
      for (uint band=0; band<bandwidth; band++)
	// there cannot be entries for offset that are too large
	if (offset(band) >= -i && offset(band) < N-i)
	  w(i)+= entry(band)*v(i+offset(band));
    return w;
  }
}

NumMatrix<double> PixelCovarianceMatrix::operator*(const NumMatrix<double>& M) const {
  // diagonal matrix
  if (bandwidth==1) {
    NumMatrix<double> M2=M;
    if (entry(0) != 1)
      M2 *= entry(0);
    return M2;
  } 
  // real banded matrix: for each column of matrix do row-wise multiplication
  // exploit limited index range for inner loop
  else {
    int N = M.getRows(), N2 = M.getColumns();
    NumMatrix<double> M2(N,N2);
    for (uint j=0; j<N2; j++) {
      boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<double> > mc((boost::numeric::ublas::matrix<double>&)M,j);
      boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<double> > m2c(M2,j);
      for (int i=0; i<N; i++)
	for (uint band=0; band<bandwidth; band++)
	  // there cannot be entries for offset that are too large
	  if (offset(band) >= -i && offset(band) < N-i)
	    m2c(i)+= entry(band)*mc(i+offset(band));
    }
    return M2;
  }
}


/// Operator for multiplication to general NumMatrix.
/// This is essentially an extension of NumMatrix to work together with
/// PixelCovarianceMatrix.
NumMatrix<double> operator*(const NumMatrix<double>& M, const PixelCovarianceMatrix& V) {
  // diagonal matrix
  if (V.getBandwidth()==1) {
    NumMatrix<double> M2=M;
    if (V(0,0) != 1)
      M2 *= V(0,0);
    return M2;
  } 
  // real banded matrix: exploit M*V = (V^T * M^T)^T
  else {
    return (V.transpose()*M.transpose()).transpose();
  }
}
