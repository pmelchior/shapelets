#include <modelfit/GalaxyModel.h>

using namespace shapelens;

GalaxyModel::~GalaxyModel() {}
void GalaxyModel::setObject(Object& obj, data_t normalization, bool add) const {
  int x,y;
  data_t x_rel, y_rel;
  for(int i=0; i < obj.size(); i++) {
    obj.grid.getCoords(i,x,y);
    x_rel = x - obj.centroid(0);  // coords must be around (0,0)
    y_rel = y - obj.centroid(1);
    if (add)
      obj(i) += getValue(x_rel,y_rel)/normalization;
    else
      obj(i) = getValue(x_rel,y_rel)/normalization;
  }
}

void GalaxyModel::setObjectSheared(Object& obj, complex<data_t> gamma, data_t normalization, bool add) const {
  NumVector<data_t> sourceCoords(2), lensCoords(2);
  NumMatrix<data_t> A(2,2);
  int x,y;
  A(0,0) = 1-real(gamma);
  A(0,1) = -imag(gamma);
  A(1,0) = -imag(gamma);
  A(1,1) = 1+real(gamma);
  for(int i=0; i < obj.size(); i++) {
    obj.grid.getCoords(i,x,y);
    lensCoords(0) = x - obj.centroid(0);
    lensCoords(1) = y - obj.centroid(1);
    sourceCoords = A*lensCoords;
    if (add)
      obj(i) += getValue(sourceCoords(0),sourceCoords(1)) * (1 - gsl_pow_2(abs(gamma))) / normalization; // shear changes size and thus also flux
    else
      obj(i) = getValue(sourceCoords(0),sourceCoords(1)) * (1 - gsl_pow_2(abs(gamma))) / normalization;
  }
}

// ##### Sersic Model ##### //
SersicModel::SersicModel(data_t n, data_t Re) : 
  n(n), Re(Re) {
  limit = 5*Re;
  b = 1.9992*n - 0.3271;
  data_t RRe1n = pow(limit/Re,1./n);
  // flux at limit
  flux_limit = exp(-b*(RRe1n -1));
  // compute total flux of model (considering the truncation at limit)
  flux = gsl_pow_2(Re)*2*M_PI*n*exp(b)/pow(b,2*n) * (gsl_sf_gamma(2*n) - gsl_sf_gamma_inc(2*n,b*RRe1n));
  // subtract level at limit such that the profile vanishes there
  flux -= M_PI*gsl_pow_2(limit)*flux_limit;
}

data_t SersicModel::getValue(data_t x, data_t y) const {
  data_t radius = sqrt(gsl_pow_2(x) + gsl_pow_2(y));
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
  // flux at limit
  flux_limit = pow(1+alpha*gsl_pow_2(limit),-beta);
  // compute total flux of model (considering the truncation at limit)
  flux = 2*M_PI*(-1 + pow(1+alpha*gsl_pow_2(limit),1-beta))/(2*alpha - 2*alpha*beta);
  // subtract level at R such that the profile vanishes there
  flux -= M_PI*gsl_pow_2(limit)*flux_limit;
}

data_t MoffatModel::getValue(data_t x, data_t y) const {
  data_t radius = sqrt(gsl_pow_2(x) + gsl_pow_2(y));
  if (radius < limit)
    return pow(1+alpha*gsl_pow_2(radius),-beta) - flux_limit;
  else
    return 0;
}

data_t MoffatModel::getFlux() const {
  return flux;
}


// ##### Interpolated Model ##### //
InterpolatedModel::InterpolatedModel(Object& obj) : obj(obj) {}

data_t InterpolatedModel::getValue(data_t x, data_t y) const {
  data_t f11, f12, f21, f22,val = 0;
  int x1,x2,y1,y2;
  x += obj.centroid(0);
  y += obj.centroid(1);
  x1 = (int) floor(x);
  x2 = x1+1;
  y1 = (int) floor(y);
  y2 = y1+1;
  if (x1 >= 0 && x2 < obj.getSize(0) && y1 >= 0 && y2 < obj.getSize(1)) {
    f11 = obj(x1,y1);
    f12 = obj(x1,y2);
    f21 = obj(x2,y1);
    f22 = obj(x2,y2);
      
    val = f11*(x2-x)*(y2-y) + f12*(x2-x)*(y-y1) + f21*(x-x1)*(y2-y) + f22*(x-x1)*(y-y1);
  }
  return val;
}

data_t InterpolatedModel::getFlux() const {
  return obj.flux;
}

// ##### Shapelet Model ##### //
ShapeletModel::ShapeletModel(ShapeletObject& sobj) : sobj(sobj) {}

data_t ShapeletModel::getValue(data_t x, data_t y) const {
  Point2D<data_t> centroid = sobj.getCentroid();
  centroid(0) += x;
  centroid(1) += y;
  return sobj.eval(centroid);
}

data_t ShapeletModel::getFlux() const {
  return sobj.getShapeletFlux();
}
