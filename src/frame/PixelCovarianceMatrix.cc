#include <frame/PixelCovarianceMatrix.h>
#include <map>
#include <vector>
#include <gsl/gsl_math.h>
#include <fstream>

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
  entry = NumVector<data_t>(bandwidth);
}

void PixelCovarianceMatrix::setBand(unsigned int band, int os, data_t val) {
  if (band<bandwidth) {
    offset(band) = os;
    entry(band) = val;
  }
}

int PixelCovarianceMatrix::getOffset(unsigned int band) const {
  if (band<bandwidth)
    return offset(band);
  else
    return 0;
}

data_t PixelCovarianceMatrix::getValue(unsigned int band) const {
  if (band<bandwidth)
    return entry(band);
  else
    return 0;
}

void PixelCovarianceMatrix::save(std::string filename) const {
  std::fstream pcm_file(filename.c_str(),std::ios::out);
  for (uint b = 0; b<bandwidth; b++)
    pcm_file << offset(b) << "\t" << entry(b) << std::endl;
  pcm_file.close();
}

void PixelCovarianceMatrix::load(std::string filename) {
  int o;
  data_t e;
  std::fstream pcm_file(filename.c_str(),std::ios::in);
  bandwidth = 0;
  offset = NumVector<int>();
  entry = NumVector<data_t>();
  while (pcm_file.good()) {
    bandwidth++;
    offset.resize(bandwidth);
    entry.resize(bandwidth);
    pcm_file >> offset(bandwidth-1) >> entry(bandwidth-1);
  }
  bandwidth--;
  offset.resize(bandwidth);
  entry.resize(bandwidth);
  pcm_file.close();
}

data_t PixelCovarianceMatrix::operator()(unsigned int i, unsigned int j) const {
  int k;
  for (uint b=0; b <bandwidth; b++) {
    k = (int)i + offset(b);
    if (k==j)
      return entry(b);
  }
  return 0;
}


void PixelCovarianceMatrix::setCovarianceMatrix(const CorrelationFunction& xi, const Grid& grid) {
  const std::map<Point2D<grid_t>, data_t>& corr = xi.getCorrelationFunction();
  const std::map<Point2D<grid_t>, data_t>& sigma = xi.getCorrelationError();
  uint N1 = grid.getSize(0);
  bandwidth = corr.size();
  offset.resize(bandwidth);
  entry.resize(bandwidth);
  uint band = 0;

  for (std::map<Point2D<grid_t>, data_t>::const_iterator iter = corr.begin(); iter !=corr.end(); iter++) {
    offset(band) = (iter->first)(0) + (iter->first)(1)*N1;
    entry(band) = iter->second;
    band++;
  }
}

NumMatrix<data_t> PixelCovarianceMatrix::getMatrix(unsigned int N) const {
  NumMatrix<data_t> M(N,N);
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
    uint N = GSL_MAX_INT(10,offset.max())*5;
    if (N%2==1) N++; // ensure that divides by 2
    NumMatrix<data_t> M = getMatrix(N).invert();
 
    int center = N/2;
    data_t max = fabs(M(center,center));
    uint inv_bandwidth = 0;
    std::vector<int> inv_offset;
    std::vector<data_t> inv_entry;
    for (int j=1; j<N-1; j++) {
      data_t val = M(N-j,j);
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
  

NumVector<data_t> PixelCovarianceMatrix::operator*(const NumVector<data_t>& v) const {
  // diagonal matrix
  if (bandwidth==1) {
    NumVector<data_t> w=v;
    if (entry(0) != 1) 
      w *= entry(0);
    return w;
  } 
  // real banded matrix: do row-wise multiplication
  // exploit limited index range for inner loop
  else {
    int N = v.size();
    NumVector<data_t> w(N);
    for (int i=0; i<N; i++)
      for (uint band=0; band<bandwidth; band++)
	// there cannot be entries for offset that are too large
	if (offset(band) >= -i && offset(band) < N-i)
	  w(i)+= entry(band)*v(i+offset(band));
    return w;
  }
}

NumMatrix<data_t> PixelCovarianceMatrix::operator*(const NumMatrix<data_t>& M) const {
  // diagonal matrix
  if (bandwidth==1) {
    NumMatrix<data_t> M2=M;
    if (entry(0) != 1)
      M2 *= entry(0);
    return M2;
  } 
  // real banded matrix: for each column of matrix do row-wise multiplication
  // exploit limited index range for inner loop
  else {
    int N = M.getRows(), N2 = M.getColumns();
    NumMatrix<data_t> M2(N,N2);
    for (uint j=0; j<N2; j++) {
      boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<data_t> > mc((boost::numeric::ublas::matrix<data_t>&)M,j);
      boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<data_t> > m2c(M2,j);
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
NumMatrix<data_t> operator*(const NumMatrix<data_t>& M, const PixelCovarianceMatrix& V) {
  // diagonal matrix
  if (V.getBandwidth()==1) {
    NumMatrix<data_t> M2=M;
    if (V(0,0) != 1)
      M2 *= V(0,0);
    return M2;
  } 
  // real banded matrix: exploit M*V = (V^T * M^T)^T
  else {
    return (V.transpose()*M.transpose()).transpose();
  }
}
