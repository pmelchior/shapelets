#include <PixelCovarianceMatrix.h>
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
  entry = NumVector<double>(bandwidth);
}

void PixelCovarianceMatrix::setBand(unsigned int band, int os, double val) {
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

double PixelCovarianceMatrix::getValue(unsigned int band) const {
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
  double e;
  std::fstream pcm_file(filename.c_str(),std::ios::in);
  bandwidth = 0;
  offset = NumVector<int>();
  entry = NumVector<double>();
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

double PixelCovarianceMatrix::operator()(unsigned int i, unsigned int j) const {
  int k;
  for (uint b=0; b <bandwidth; b++) {
    k = (int)i + offset(b);
    if (k==j)
      return entry(b);
  }
  return 0;
}


void PixelCovarianceMatrix::setCovarianceMatrix(const CorrelationFunction& xi, const Grid& grid, History& hist) {
  hist.append("# Reading correlation function:\n");
  const NumVector<double>& corr = xi.getCorrelationFunction();
  const NumVector<double>& distance = xi.getDistances();
  const NumVector<double>& sigma = xi.getCorrelationError();
  if (corr.size() == 0) {
    std::cout << "PixelCovarianceMatrix: correlation function has size 0" << std::endl;
    std::terminate();
  }
  // print correlation function in a nice way
  hist.append("# dist\tcorr\t\t +- sigma\n");
  char diststr[10];
  char corrstr[20];
  for (int i=0; i < corr.size(); i++) {
    sprintf(diststr,"%1.2f",distance(i)); 
    sprintf(corrstr,"%1.2e\t +- %1.2e",corr(i),sigma(i));
    hist.append("# "+std::string(diststr)+"\t"+std::string(corrstr)+"\n");
  }

  // the bandwidth depends on the values of the correlation function
  // for distances larger than 0.
  // use 2% of the maximum (distance 0) value as threshold, or the information content 
  // of the correlation function to decide what is necessary for proper description
  uint N1 = grid.getSize(0);
  // effectively zero pixel correlation
  if (corr.size() == 1 || corr(1) < 0.02*corr(0)) {
    bandwidth = 1;
    hist.append("# Restricting covariance matrix to 1 band: corr(r>0) = 0\n");
    setBandwidth(bandwidth);
    offset(0) = 0;
    entry(0) = corr(0);
  }
  // include all pixels with distance 1
  else if (corr.size() >= 3 && corr(2) < 0.02*corr(0)) {
    bandwidth = 5;
    hist.append("# Restricting covariance matrix to 5 band: corr(r>1) = 0\n");
    setBandwidth(bandwidth);
    offset(0) = -N1;
    offset(1) = -1;
    offset(2) = 0;
    offset(3) = 1;
    offset(4) = N1;
    entry(2) = corr(0);
    entry(0) = entry(1) = entry(3) = entry(4) = corr(1);
  }
  // include now also top-right/top-left/bottom-right/bottom-left neighbor: 3x3 box
  else if (corr.size() == 3 || corr(3) < 0.02*corr(0)){
    bandwidth = 9;
    hist.append("# Restricting covariance matrix to 9 bands: corr(r>=2) = 0\n");
    setBandwidth(bandwidth);
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
  // distance = 2 included
  else if (corr.size() >= 6 && corr(4) < 0.02*corr(0)) {
    bandwidth = 13;
    hist.append("# Restricting covariance matrix to 13 bands: corr(r>2) = 0\n");
    setBandwidth(bandwidth);
    offset(0) = -2*N1;
    offset(1) = -N1-1;
    offset(2) = -N1;
    offset(3) = -N1+1;
    offset(4) = -2;
    offset(5) = -1;
    offset(6) = 0;
    offset(7) = 1;
    offset(8) = 2;
    offset(9) = N1-1;
    offset(10) = N1;
    offset(11) = N1+1;
    offset(12) = 2*N1;
    entry(0) = entry(4) = entry(8) = entry(12) = corr(3);
    entry(1) = entry(3) = entry(9) = entry(11) = corr(2);
    entry(2) = entry(5) = entry(7) = entry(10) = corr(1);
    entry(6) = corr(0);
  }
  // distance = sqrt(5) included
  else if (corr.size() >= 6 && corr(5) < 0.02*corr(0)) {
    bandwidth = 21;
    hist.append("# Restricting covariance matrix to 21 bands: corr(r>sqrt(5)) = 0\n");
    setBandwidth(bandwidth);
    offset(0) = -2*N1-1;
    offset(1) = -2*N1;
    offset(2) = -2*N1+1;
    offset(3) = -N1-2;
    offset(4) = -N1-1;
    offset(5) = -N1;
    offset(6) = -N1+1;
    offset(7) = -N1+2;
    offset(8) = -2;
    offset(9) = -1;
    offset(10) = 0;
    offset(11) = 1;
    offset(12) = 2;
    offset(13) = N1-2;
    offset(14) = N1-1;
    offset(15) = N1;
    offset(16) = N1+1;
    offset(17) = N1+2;
    offset(18) = 2*N1-1;
    offset(19) = 2*N1;
    offset(20) = 2*N1+1;
    entry(0) = entry(2) = entry(3) = entry(7) = entry(13) = entry(17) = entry(18) = entry(20) = corr(4);
    entry(1) = entry(8) = entry(12) = entry(19) = corr(3);
    entry(4) = entry(6) = entry(14) = entry(16) = corr(2);
    entry(5) = entry(9) = entry(11) = entry(15) = corr(1);
    entry(10) = corr(0);
  }
  // include all pixels in a box of 5x5
  else {
    bandwidth = 25;
    hist.append("# Restricting covariance matrix to 25 bands: corr(r>=3)) = 0\n");
    setBandwidth(bandwidth);
    offset(0) = -2*N1-2;
    offset(1) = -2*N1-1;
    offset(2) = -2*N1;
    offset(3) = -2*N1+1;
    offset(4) = -2*N1+2;
    offset(5) = -N1-2;
    offset(6) = -N1-1;
    offset(7) = -N1;
    offset(8) = -N1+1;
    offset(9) = -N1+2;
    offset(10) = -2;
    offset(11) = -1;
    offset(12) = 0;
    offset(13) = 1;
    offset(14) = 2;
    offset(15) = N1-2;
    offset(16) = N1-1;
    offset(17) = N1;
    offset(18) = N1+1;
    offset(19) = N1+2;
    offset(20) = 2*N1-2;
    offset(21) = 2*N1-1;
    offset(22) = 2*N1;
    offset(23) = 2*N1+1;
    offset(24) = 2*N1+2;
    entry(0) = entry(4) = entry(20) = entry(24) = corr(5);
    entry(1) = entry(3) = entry(5) = entry(9) = entry(15) = entry(19) = entry(21) = entry(23) = corr(4);
    entry(2) = entry(10) = entry(14) = entry(22) = corr(3);
    entry(6) = entry(8) = entry(16) = entry(18) = corr(2);
    entry(7) = entry(11) = entry(13) = entry(17) = corr(1);
    entry(12) = corr(0);
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
    uint N = GSL_MAX_INT(10,offset.max())*5;
    if (N%2==1) N++; // ensure that divides by 2
    NumMatrix<double> M = getMatrix(N).svd_invert();
    writeFITSFile("V.fits",getMatrix(N));
    addFITSExtension("V.fits","SVD_INV",M);

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
