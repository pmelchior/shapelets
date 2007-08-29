#include <PolarTransformation.h>
#include <MatrixManipulations.h>
#include <NumVector.h>

typedef complex<double> Complex;
const Complex I = Complex(0,1);

PolarTransformation::PolarTransformation() {
}

PolarTransformation::PolarTransformation(unsigned int innmax) {
  nmax = innmax;
  nVector = IndexVector(nmax);
  int nCoeffs = nVector.getNCoeffs();
  c2p = NumMatrix<Complex>(nCoeffs,nCoeffs);
  p2c = NumMatrix<Complex>(nCoeffs,nCoeffs);
  buildTransformationMatrix();
}

unsigned int PolarTransformation::getOrder() {
  return nmax;
}

void PolarTransformation::setOrder (unsigned int innmax) {
  if (innmax != nmax) {
    nmax = innmax;
    nVector = IndexVector(nmax);
    int nCoeffs = nVector.getNCoeffs();
    c2p = NumMatrix<Complex>(nCoeffs,nCoeffs);
    buildTransformationMatrix();
  }
}

void PolarTransformation::getPolarCoeffs(const NumMatrix<double>& cartesianCoeffs, NumMatrix<Complex>& polarCoeffs) {
  NumVector<Complex> cartesianVector, polarVector;
  // convert diagonalized form of matrix into vector, implicit typecast to Complex
  matrixMapping(cartesianCoeffs,cartesianVector,0,nVector);
  // the actual transformation in reduced coeff space
  polarVector = c2p*cartesianVector;
  // reconstruct the matrix from the vector, in polar matrix form
  vectorMapping(polarVector,polarCoeffs,nVector);
}
void PolarTransformation::getCartesianCoeffs(const NumMatrix<Complex>& polarCoeffs, NumMatrix<double>& cartesianCoeffs){
  NumVector<Complex> polarVector, cartesianVector;
  matrixMapping(polarCoeffs,polarVector,1,nVector);
  cartesianVector = p2c*polarVector;
  // this time reconstruct into cartesian form, aka upper left triangular matrix form
  // explicit typecast from Complex to double is performed
  vectorMapping(cartesianVector,cartesianCoeffs,nVector);
}


// see Paper III, eq. 14
void PolarTransformation::buildTransformationMatrix() {
  Complex prefactor, factorials, doublesum;
  unsigned int n_;
  int m_;
  int nCoeffs = nVector.getNCoeffs();
  for(int i = 0; i < nCoeffs; i++) {     // polar vector index
    for (int j = 0; j < nCoeffs; j++) {  // cartesian vector index
      // determine n,m,n1,n2 from i,j
      int nr = nVector.getN1(i), nl = nVector.getN2(i);
      int n1 = nVector.getN1(j), n2 = nVector.getN2(j);
      // the prefactor
      prefactor = pow(2,-(double)(nr+nl)/2)*pow(I,nr-nl);
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
	      doublesum += pow(I,- m_) * gsl_sf_fact(nr) * gsl_sf_fact(nl) /
		(gsl_sf_fact(nr_)*gsl_sf_fact(nr - nr_) * 
		 gsl_sf_fact(nl_)*gsl_sf_fact(nl - nl_));
	    }
	  }
	}
      } else factorials = 0;
      c2p(i,j) = prefactor * factorials * doublesum ;
    }
  }
  p2c = c2p.invert();
}
      
