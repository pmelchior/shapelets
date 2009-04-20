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

void SourceModel::setEllipticalSupport(data_t radius, const complex<data_t>& eps) {
  // compute orientation angle
  data_t theta;
  if (real(eps)!=0)
    theta = 0.5*atan(imag(eps)/real(eps));
  else
    theta = 0;
  // theta can not distinguish between vertical and horizontal orientation:
  // theta = 0 in both cases
  // since we want to have angle to x-Axis, map theta onto unbroken range of 180 deg
  if (real(eps) < 0)
    theta += M_PI_2;
  // compute size of semi-major and semi-minor axis (axis-parallel system)
  data_t a = (1 + fabs(real(eps)))*radius;
  data_t b = (1 - fabs(real(eps)))*radius;
  // compute curve parameter t which maximizes x or y
  data_t tx = atan(-b*tan(theta)/a);
  data_t ty = atan(b/(tan(theta)*a));
  // insert in parametric equation for rotated ellipse
  data_t max_x = fabs(a*cos(tx)*cos(theta) - b*sin(tx)*sin(theta));
  data_t max_y = fabs(a*cos(ty)*sin(theta) + b*sin(ty)*cos(theta));
  support.tr(0) = max_x;
  support.tr(1) = max_y;
  support.ll(0) = - support.tr(0);
  support.ll(1) = - support.tr(1);
  support.ll += centroid;
  support.tr += centroid;
}

// ##### Sersic Model ##### //
SersicModel::SersicModel(data_t n, data_t Re, complex<data_t> eps, const Point2D<data_t>& centroid) : 
  n(n), Re(Re), eps(eps) {
  limit = 5*Re;
  shear_norm = 1 - gsl_pow_2(abs(eps));
  SourceModel::centroid = centroid;
  SourceModel::setEllipticalSupport(limit,eps);
  b = 1.9992*n - 0.3271;
  data_t RRe1n = pow(limit/Re,1./n);
  // flux at limit
  flux_limit = exp(-b*(RRe1n -1));
  // compute total flux of model (considering the truncation at limit)
  flux = (gsl_pow_2(Re)*2*M_PI*n*exp(b)/pow(b,2*n) * (gsl_sf_gamma(2*n) - gsl_sf_gamma_inc(2*n,b*RRe1n)));
  // subtract level at limit such that the profile vanishes there
  flux -= M_PI*gsl_pow_2(limit)*flux_limit;
  // correct for shearing
  flux /= shear_norm;
}

data_t SersicModel::getValue(const Point2D<data_t>& P) const {
  data_t x = P(0) - centroid(0);
  data_t y = P(1) - centroid(1);
  data_t x_ = (1-real(eps))*x - imag(eps)*y;
  data_t y_ = -imag(eps)*x + (1+real(eps))*y;
  
  data_t radius = sqrt(x_*x_ + y_*y_)/shear_norm; // shear changes size
  if (radius < limit)
    return exp(-b*(pow(radius/Re,1./n) -1)) - flux_limit;
  else
    return 0;
}

data_t SersicModel::getFlux() const {
  return flux;
}

// ##### Moffat Model ##### //
MoffatModel::MoffatModel(data_t beta, data_t FWHM, complex<data_t> eps, const Point2D<data_t>& centroid) :
  beta(beta), eps(eps) {
  SourceModel::centroid = centroid;
  alpha = (pow(2.,1./beta)-1)/gsl_pow_2(FWHM/2);
  limit = 5*FWHM;
  shear_norm = 1 - gsl_pow_2(abs(eps));
  // define area of support
  SourceModel::centroid = centroid;
  SourceModel::setEllipticalSupport(limit,eps);
  // flux at limit
  flux_limit = pow(1+alpha*gsl_pow_2(limit),-beta);
  // compute total flux of model (considering the truncation at limit)
  flux = 2*M_PI*(-1 + pow(1+alpha*gsl_pow_2(limit),1-beta))/(2*alpha - 2*alpha*beta);
  // subtract level at R such that the profile vanishes there
  flux -= M_PI*gsl_pow_2(limit)*flux_limit;
  // correct for shearing
  flux /= shear_norm;
}

data_t MoffatModel::getValue(const Point2D<data_t>& P) const {
  data_t x = P(0) - centroid(0);
  data_t y = P(1) - centroid(1);
  data_t x_ = (1-real(eps))*x - imag(eps)*y;
  data_t y_ = -imag(eps)*x + (1+real(eps))*y;
  
  data_t radius = sqrt(x_*x_ + y_*y_)/shear_norm;
  if (radius < limit)
    return pow(1+alpha*gsl_pow_2(radius),-beta) - flux_limit;
  else
    return 0;
}

data_t MoffatModel::getFlux() const {
  return flux;
}


// ##### Interpolated Model ##### //
InterpolatedModel::InterpolatedModel(const Object& obj, int order) : 
  obj(obj), order(order) {
  // define area of support
  support.ll = Point2D<data_t>(obj.grid.getStartPosition(0),obj.grid.getStartPosition(1));
  support.tr = Point2D<data_t>(obj.grid.getStopPosition(0),obj.grid.getStopPosition(1));
}

data_t InterpolatedModel::getValue(const Point2D<data_t>& P) const {
  switch (order) {
  case 1: // simple bi-linear interpolation
    return obj.interpolate(P(0)-centroid(0),P(1)-centroid(1));
  case -3: // bi-cubic interpolation
    return Interpolation::bicubic(obj,P(0)-centroid(0),P(1)-centroid(1));
  default: // nth-order polynomial interpolation
    return Interpolation::polynomial(obj,P(0)-centroid(0),P(1)-centroid(1),order);
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
  return sobj.eval(centroid);
}

data_t ShapeletModel::getFlux() const {
  return sobj.getShapeletFlux();
}
