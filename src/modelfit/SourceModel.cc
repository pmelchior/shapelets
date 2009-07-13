#include "../../include/modelfit/SourceModel.h"
#include "../../include/utils/Interpolation.h"

using namespace shapelens;

SourceModel::~SourceModel() {}

const Rectangle<data_t>& SourceModel::getSupport() const {
  return support;
}

const Point<data_t>& SourceModel::getCentroid() const {
  return centroid;
}

unsigned long SourceModel::getID() const {
  return id;
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

// ##### SourceModelList ##### //
Catalog SourceModelList::getCatalog() const {
  Catalog cat;
  CatObject co;
  co.FLAGS = 0;
  for (unsigned long i=0; i < SourceModelList::size(); i++) {
    const SourceModel& sm = *SourceModelList::operator[](i);
    const Rectangle<data_t>& support = sm.getSupport();
    const Point<data_t>& centroid = sm.getCentroid();
    co.XMIN = support.ll(0);
    co.YMIN = support.ll(1);
    co.XMAX = support.tr(0);
    co.YMAX = support.tr(1);
    co.XCENTROID = centroid(0);
    co.YCENTROID = centroid(1);
    co.FLUX = sm.getFlux();
    co.CLASSIFIER = sm.getModelType();
    co.PARENT = sm.getID();
    cat[i] = co;
  }
  return cat;
}


// ##### Sersic Model ##### //
SersicModel::SersicModel(data_t n, data_t Re, data_t flux_eff, complex<data_t> eps, const Point<data_t>& centroid, unsigned long id) : 
  n(n), Re(Re), eps(eps) {
  limit = 5*Re;
  shear_norm = 1 - gsl_pow_2(abs(eps));
  SourceModel::centroid = centroid;
  SourceModel::setEllipticalSupport(limit,eps);
  SourceModel::id = id;
  b = 1.9992*n - 0.3271;
  data_t RRe1n = pow(limit/Re,1./n);
  // flux at limit
  flux_limit = exp(-b*(RRe1n -1));
  // compute total flux of model (considering the truncation at limit)
  flux = (gsl_pow_2(Re)*2*M_PI*n*exp(b)/pow(b,2*n) * (gsl_sf_gamma(2*n) - gsl_sf_gamma_inc(2*n,b*RRe1n)));
  // subtract level at limit such that the profile vanishes there
  flux -= M_PI*gsl_pow_2(limit)*flux_limit;
  // correct for shearing
  flux *= shear_norm;
  // compute rescaling factor for flux
  flux_scale = flux_eff/flux;
}

data_t SersicModel::getValue(const Point<data_t>& P) const {
  data_t x = P(0) - centroid(0);
  data_t y = P(1) - centroid(1);
  data_t x_ = (1-real(eps))*x - imag(eps)*y;
  data_t y_ = -imag(eps)*x + (1+real(eps))*y;
  
  data_t radius = sqrt(x_*x_ + y_*y_)/shear_norm; // shear changes size
  if (radius < limit)
    return flux_scale*(exp(-b*(pow(radius/Re,1./n) -1)) - flux_limit);
  else
    return 0;
}

data_t SersicModel::getFlux() const {
  return flux_scale*flux;
}

char SersicModel::getModelType() const {
  return 0;
}

// ##### Moffat Model ##### //
MoffatModel::MoffatModel(data_t beta, data_t FWHM, data_t flux_eff, complex<data_t> eps, const Point<data_t>& centroid, unsigned long id) :
  beta(beta), eps(eps) {
  SourceModel::centroid = centroid;
  SourceModel::id = id;
  alpha = (pow(2.,1./beta)-1)/gsl_pow_2(FWHM/2);
  limit = 5*FWHM;
  shear_norm = 1 - gsl_pow_2(abs(eps));
  // define area of support
  SourceModel::setEllipticalSupport(limit,eps);
  // flux at limit
  flux_limit = pow(1+alpha*gsl_pow_2(limit),-beta);
  // compute total flux of model (considering the truncation at limit)
  flux = 2*M_PI*(-1 + pow(1+alpha*gsl_pow_2(limit),1-beta))/(2*alpha - 2*alpha*beta);
  // subtract level at R such that the profile vanishes there
  flux -= M_PI*gsl_pow_2(limit)*flux_limit;
  // correct for shearing
  flux *= shear_norm;
  // compute rescaling factor for flux
  flux_scale = flux_eff/flux;
}

data_t MoffatModel::getValue(const Point<data_t>& P) const {
  data_t x = P(0) - centroid(0);
  data_t y = P(1) - centroid(1);
  data_t x_ = (1-real(eps))*x - imag(eps)*y;
  data_t y_ = -imag(eps)*x + (1+real(eps))*y;
  
  data_t radius = sqrt(x_*x_ + y_*y_)/shear_norm;
  if (radius < limit)
    return flux_scale*(pow(1+alpha*gsl_pow_2(radius),-beta) - flux_limit);
  else
    return 0;
}

data_t MoffatModel::getFlux() const {
  return flux_scale*flux;
}

char MoffatModel::getModelType() const {
  return 3;
}

// ##### Interpolated Model ##### //
InterpolatedModel::InterpolatedModel(const Image<data_t>& im, data_t flux, const Point<data_t>& reference, int order, unsigned long id) : 
  im(im), order(order),flux(flux), reference(reference) {
  // define area of support
  SourceModel::support = im.grid.getBoundingBox();
  // no need to defined centroid, as our reference is left-lower corner
  SourceModel::support.ll -= reference;
  SourceModel::support.tr -= reference;
  SourceModel::id = id;
  // compute total flux from pixel values
  data_t f = 0;
  for (unsigned long i =0; i < im.size(); i++)
    f += im(i);
  // compute rescaling factor for flux
  flux_scale = flux/f;
}

data_t InterpolatedModel::getValue(const Point<data_t>& P) const {
  Point<data_t> P_ = P;
  P_ -= reference;
  switch (order) {
  case 1: // simple bi-linear interpolation
    return flux_scale*im.interpolate(P_);
  case -3: // bi-cubic interpolation
    return flux_scale*Interpolation::bicubic(im,P_);
  default: // nth-order polynomial interpolation
    return flux_scale*Interpolation::polynomial(im,P_,order);
  }
}

data_t InterpolatedModel::getFlux() const {
  return flux;
}

char InterpolatedModel::getModelType() const {
  return 2;
}

// ##### Shapelet Model ##### //
ShapeletModel::ShapeletModel(const ShapeletObject& sobj, data_t flux, const Point<data_t>& centroid) : 
  sobj(sobj), scentroid(sobj.getCentroid()) {
  SourceModel::centroid = centroid;
  SourceModel::id = sobj.getObjectID();
  // define area of support
  // adjust for new location of centroid
  const Grid& grid = sobj.getGrid();
  support.ll = Point<data_t>(grid.getStartPosition(0)-scentroid(0)+centroid(0),grid.getStartPosition(1)-scentroid(1)+centroid(1));
  support.tr = Point<data_t>(grid.getStopPosition(0)-scentroid(0)+centroid(0)-1,grid.getStopPosition(1)-scentroid(1)+centroid(1)-1);
  flux_scale = flux/sobj.getShapeletFlux();
}

data_t ShapeletModel::getValue(const Point<data_t>& P) const {
  Point<data_t> P_ = P;
  P_ -= centroid;
  P_ += scentroid;
  return flux_scale*sobj.eval(P_);
}

data_t ShapeletModel::getFlux() const {
  return flux_scale*sobj.getShapeletFlux();
}

char ShapeletModel::getModelType() const {
  return 1;
}
