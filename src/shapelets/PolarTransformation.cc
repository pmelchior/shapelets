#include <shapelets/PolarTransformation.h>
#include <shapelets/CoefficientVector.h>

typedef complex<data_t> Complex;
const Complex I = Complex(0,1);

PolarTransformation::PolarTransformation() {
}

void PolarTransformation::getPolarCoeffs(const CoefficientVector<data_t>& cartesianCoeffs, CoefficientVector<Complex>& polarCoeffs, NumMatrix<data_t>* cov, NumMatrix<complex<data_t> >* polarCov) {
  // set nmax of polarCoeffs to nmax of cartesianCoeffs
  polarCoeffs.setNMax(cartesianCoeffs.getNMax());
  // set up c2p; automatically adjusts to order of input coeffs
  buildTransformationMatrix(cartesianCoeffs.getIndexVector(),polarCoeffs.getIndexVector());
  // intermediate step: as c2p is complex, we need complex factors
  NumVector<Complex> cartV = cartesianCoeffs.getNumVector();
  // the actual transformation in shapelet space
  polarCoeffs = c2p*cartV;
  // apply transformation to covariance matrices
  if (cov != NULL && polarCov != NULL && cov->getColumns() == cartesianCoeffs.getNCoeffs()) {
    NumMatrix<Complex> covM = *cov;
    *polarCov = (c2p*covM)*c2p.transpose();
  }
}

void PolarTransformation::getCartesianCoeffs(const CoefficientVector<Complex>& polarCoeffs, CoefficientVector<data_t>& cartesianCoeffs, NumMatrix<complex<data_t> >* polarCov, NumMatrix<data_t>* cov){
  // set nmax of polarCoeffs to nmax of cartesianCoeffs
  cartesianCoeffs.setNMax(polarCoeffs.getNMax());
  // set up p2c; automatically adjusts to order of input coeffs
  buildTransformationMatrix(cartesianCoeffs.getIndexVector(),polarCoeffs.getIndexVector());
  // intermediate step necessary since p2c is complex
  NumVector<Complex> cartV = p2c*polarCoeffs;
  cartesianCoeffs.setNMax(polarCoeffs.getNMax());
  // imaginary part of cartesianVector should be small
  for (unsigned int i=0; i<cartesianCoeffs.size(); i++)
    cartesianCoeffs(i) = std::real(cartV(i));
  if (cov != NULL && polarCov != NULL && polarCov->getColumns() == polarCoeffs.getNCoeffs()) {
    NumMatrix<Complex> covM = (p2c*(*polarCov))*p2c.transpose();
    for (unsigned int i=0; i<covM.getRows(); i++)
      for (unsigned int j=0; j<covM.getColumns(); j++)
	(*cov)(i,j) = std::real(covM(i,j));
  }
}


// see Paper III, eq. 14
void PolarTransformation::buildTransformationMatrix(const IndexVector& nVector, const IndexVector& pnVector) {
  int nCoeffs = nVector.getNCoeffs();
  // only do something if orders do not match
  if (c2p.getRows() != nCoeffs) {
    c2p.resize(nCoeffs,nCoeffs);
    c2p.clear();
 
    Complex prefactor, factorials, doublesum;
    int n,n_,m,m_,n1,n2,nr,nl;
    for(int i = 0; i < nCoeffs; i++) {     // polar vector index
      n = pnVector.getState1(i);
      m = pnVector.getState2(i);
      nr = (n+m)/2;
      nl = (n-m)/2;
      for (int j = 0; j < nCoeffs; j++) {  // cartesian vector index
	n1 = nVector.getState1(j);
	n2 = nVector.getState2(j);
	// the prefactor
	prefactor = data_t(pow(2,-data_t(nr+nl)/2))*pow(I,nr-nl);
	// only orders with n1+n2 <= nmax are evaluated
	// since the cartesian orders are symmetric, n1 and n2 are bounded by nmax
	//if (n1+n2 <= nmax) {
	// the first kronecker
	if ((n1+n2) == (nr+nl)) {
	  factorials = sqrt(gsl_sf_fact(n1) * gsl_sf_fact(n2) / 
			    (gsl_sf_fact(nr)*gsl_sf_fact(nl)) );
	  // the double sum
	  doublesum = 0;
	  for (int nr_ = 0; nr_ <= nr; nr_++) {
	    for (int nl_ = 0; nl_ <= nl; nl_++) {
	      n_ = nr_ + nl_;
	      m_ = nr_ - nl_;
	      // the second kronecker
	      if (n_ == n1) {
		doublesum += data_t(gsl_sf_fact(nr) * gsl_sf_fact(nl) /
				    (gsl_sf_fact(nr_)*gsl_sf_fact(nr - nr_) * 
				     gsl_sf_fact(nl_)*gsl_sf_fact(nl - nl_))) * pow(I,- m_);
	      }
	    }
	  }
	} else factorials = 0;
	c2p(i,j) = prefactor * factorials * doublesum ;
      }
    }
    p2c = c2p.invert();
  }
}
      
