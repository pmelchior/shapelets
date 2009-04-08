#include <modelfit/SourceModel.h>
#include <utils/Interpolation.h>

using namespace shapelens;

SourceModel::~SourceModel() {}

const Rectangle<data_t>& SourceModel::getSupport() const {
  return support;
}

void SourceModel::setObject(Object& obj, data_t normalization, bool add) const {
  int x,y;
  Point2D<data_t> P;
  for(int i=0; i < obj.size(); i++) {
    obj.grid.getCoords(i,x,y);
    P(0) = x - obj.centroid(0);  // coords must be around (0,0)
    P(1) = y - obj.centroid(1);
    if (add)
      obj(i) += getValue(P)/normalization;
    else
      obj(i) = getValue(P)/normalization;
  }
}

void SourceModel::setObjectSheared(Object& obj, complex<data_t> gamma, data_t normalization, bool add) const {
  NumVector<data_t> sourceCoords(2), lensCoords(2);
  NumMatrix<data_t> A(2,2);
  int x,y;
  Point2D<data_t> P;
  A(0,0) = 1-real(gamma);
  A(0,1) = -imag(gamma);
  A(1,0) = -imag(gamma);
  A(1,1) = 1+real(gamma);
  for(int i=0; i < obj.size(); i++) {
    obj.grid.getCoords(i,x,y);
    lensCoords(0) = x - obj.centroid(0);
    lensCoords(1) = y - obj.centroid(1);
    sourceCoords = A*lensCoords;
    P(0) = sourceCoords(0);
    P(1) = sourceCoords(1);
    if (add)
      obj(i) += getValue(P) * (1 - gsl_pow_2(abs(gamma))) / normalization; // shear changes size and thus also flux
    else
      obj(i) = getValue(P) * (1 - gsl_pow_2(abs(gamma))) / normalization;
  }
}

// ##### Sersic Model ##### //
SersicModel::SersicModel(data_t n, data_t Re) : 
  n(n), Re(Re) {
  limit = 5*Re;
  // define area of support
  support.ll(0) = support.ll(1) = -limit;
  support.tr(0) = support.tr(1) = limit;
  b = 1.9992*n - 0.3271;
  data_t RRe1n = pow(limit/Re,1./n);
  // flux at limit
  flux_limit = exp(-b*(RRe1n -1));
  // compute total flux of model (considering the truncation at limit)
  flux = gsl_pow_2(Re)*2*M_PI*n*exp(b)/pow(b,2*n) * (gsl_sf_gamma(2*n) - gsl_sf_gamma_inc(2*n,b*RRe1n));
  // subtract level at limit such that the profile vanishes there
  flux -= M_PI*gsl_pow_2(limit)*flux_limit;
}

data_t SersicModel::getValue(const Point2D<data_t>& P) const {
  data_t radius = sqrt(gsl_pow_2(P(0)) + gsl_pow_2(P(1)));
  if (radius < limit)
    return exp(-b*(pow(radius/Re,1./n) -1)) - flux_limit;
  else
    return 0;
}

data_t SersicModel::getFlux() const {
  return flux;
}

// ##### Moffat Model ##### //
MoffatModel::MoffatModel(data_t beta, data_t FWHM) :
  beta(beta) {
  alpha = (pow(2.,1./beta)-1)/gsl_pow_2(FWHM/2);
  limit = 5*FWHM;
  // define area of support
  support.ll(0) = support.ll(1) = -limit;
  support.tr(0) = support.tr(1) = limit;
  // flux at limit
  flux_limit = pow(1+alpha*gsl_pow_2(limit),-beta);
  // compute total flux of model (considering the truncation at limit)
  flux = 2*M_PI*(-1 + pow(1+alpha*gsl_pow_2(limit),1-beta))/(2*alpha - 2*alpha*beta);
  // subtract level at R such that the profile vanishes there
  flux -= M_PI*gsl_pow_2(limit)*flux_limit;
}

data_t MoffatModel::getValue(const Point2D<data_t>& P) const {
  data_t radius = sqrt(gsl_pow_2(P(0)) + gsl_pow_2(P(1)));
  if (radius < limit)
    return pow(1+alpha*gsl_pow_2(radius),-beta) - flux_limit;
  else
    return 0;
}

data_t MoffatModel::getFlux() const {
  return flux;
}


// ##### Interpolated Model ##### //
InterpolatedModel::InterpolatedModel(Object& obj, int order) : obj(obj), order(order) {
  // define area of support
  support.ll = Point2D<data_t>(obj.grid.getStartPosition(0),obj.grid.getStartPosition(1));
  support.tr = Point2D<data_t>(obj.grid.getStopPosition(0),obj.grid.getStopPosition(1));
}

data_t InterpolatedModel::getValue(const Point2D<data_t>& P) const {
  switch (order) {
  case 1: // simple bi-linear interpolation
    return obj.interpolate(P(0)+obj.centroid(0),P(1)+obj.centroid(1));
  case -3: // bi-cubic interpolation
    return Interpolation::bicubic(obj,P(0)+obj.centroid(0),P(1)+obj.centroid(1));
  default: // nth-order polynomial interpolation
    return Interpolation::polynomial(obj,P(0)+obj.centroid(0),P(1)+obj.centroid(1),order);
  }
}

data_t InterpolatedModel::getFlux() const {
  return obj.flux;
}

// ##### Shapelet Model ##### //
ShapeletModel::ShapeletModel(const ShapeletObject& sobj) : sobj(sobj) {
  // define area of support
  const Grid& grid = sobj.getGrid();
  support.ll = Point2D<data_t>(grid.getStartPosition(0),grid.getStartPosition(1));
  support.tr = Point2D<data_t>(grid.getStopPosition(0),grid.getStopPosition(1));
}

data_t ShapeletModel::getValue(const Point2D<data_t>& P) const {
  Point2D<data_t> centroid = sobj.getCentroid();
  centroid += P;
  return sobj.eval(centroid);
}

data_t ShapeletModel::getFlux() const {
  return sobj.getShapeletFlux();
}
