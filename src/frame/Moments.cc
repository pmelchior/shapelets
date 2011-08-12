#include "../../include/frame/Moments.h"
#include "../../include/utils/MathHelper.h"

namespace shapelens {

  int SUBPIXEL = 2;
  data_t offset = 1./SUBPIXEL;

  Moment0::Moment0() {
    M = 0;
  }

  Moment0::Moment0(const Object& obj, const WeightFunction& w) {
    const Grid& grid = obj.grid;
    M = 0;
    data_t w_;
    for (long i=0; i< grid.size(); i++) {
      //Point<data_t> P = obj.grid(i), P_;
      //for (int n1 = 0; n1 < SUBPIXEL; n1++) {
      // P_(0) = P(0) + (0.5+n1)*offset - 0.5;
      //       for (int n2 = 0; n2 < SUBPIXEL; n2++) {
      // 	P_(1) = P(1) + (0.5+n2)*offset - 0.5;
      // 	w_ = w(P_)/(SUBPIXEL*SUBPIXEL);
      // 	if (obj.weight.size() != 0)
      // 	  w_ *= obj.weight(i);
      // 	M += w_ * obj.interpolate(P_);//Interpolation::bicubic(obj,P_);
      //       }
      //     }
    
      w_ = w(grid(i));
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);
      M += obj(i) * w_;
    }
  }

  data_t& Moment0::operator()(bool i) {
    return M;
  }
  const data_t& Moment0::operator()(bool i) const {
    return M;
  }

  Moment1::Moment1() {
    M[0] = M[1] = 0;
  }
  Moment1::Moment1(const Object& obj, const WeightFunction& w) {
    const Grid& grid = obj.grid;
    M[0] = M[1] = 0;
    data_t w_, diff_x, diff_y;
    Point<data_t> centroid_wcs = obj.centroid;
    grid.getWCS().transform(centroid_wcs);

    for (long i=0; i< grid.size(); i++) {
      w_ = w(grid(i));
      diff_x = obj.grid(i,0)-centroid_wcs(0);
      diff_y = obj.grid(i,1)-centroid_wcs(1);
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);
      M[0] += diff_x * obj(i) * w_;
      M[1] += diff_y * obj(i) * w_;
    }
  }
  data_t& Moment1::operator()(bool i) {
    return M[i];
  }
  const data_t& Moment1::operator()(bool i) const {
    return M[i];
  }



  Moment2::Moment2() {
    M[0] = M[1] = M[2] = 0;
  }
  Moment2::Moment2(const Object& obj, const WeightFunction& w) {
    const Grid& grid = obj.grid;
    M[0] = M[1] = M[2] = 0;
    data_t w_, diff_x, diff_y, val;
    Point<data_t> centroid_wcs = obj.centroid;
    grid.getWCS().transform(centroid_wcs);

    for (long i=0; i< grid.size(); i++) {
      val = obj(i);
      w_ = w(obj.grid(i));
      diff_x = obj.grid(i,0)-centroid_wcs(0);
      diff_y = obj.grid(i,1)-centroid_wcs(1);
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);

      M[0] += diff_x * diff_x * obj(i) * w_;
      M[1] += diff_x * diff_y * obj(i) * w_;
      M[2] += diff_y * diff_y * obj(i) * w_;
    }
  }
  data_t& Moment2::operator()(bool i, bool j) {
    return M[i+j];
  }
  const data_t& Moment2::operator()(bool i, bool j) const {
    return M[i+j];
  }



  Moment3::Moment3() {
    M[0] = M[1] = M[2] = M[3] = 0;
  }
  Moment3::Moment3(const Object& obj, const WeightFunction& w) {
    const Grid& grid = obj.grid;
    M[0] = M[1] = M[2] = M[3] = 0;
    data_t w_, diff_x, diff_y, val;
    Point<data_t> centroid_wcs = obj.centroid;
    grid.getWCS().transform(centroid_wcs);

    for (long i=0; i< grid.size(); i++) {
      val = obj(i);
      w_ = w(obj.grid(i));
      diff_x = obj.grid(i,0)-centroid_wcs(0);
      diff_y = obj.grid(i,1)-centroid_wcs(1);
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);

      M[0] += pow_int(diff_x,3) * obj(i) * w_;
      M[1] += diff_x * diff_x * diff_y * obj(i) * w_;
      M[2] += diff_x * diff_y * diff_y * obj(i) * w_;
      M[3] += pow_int(diff_y,3) * obj(i) * w_;
    }
  }
  data_t& Moment3::operator()(bool i, bool j, bool k) {
    return M[i+j+k];
  }
  const data_t& Moment3::operator()(bool i, bool j, bool k) const {
    return M[i+j+k];
  }



  Moment4::Moment4() {
    M[0] = M[1] = M[2] = M[3] = M[4] = 0;
  }
  Moment4::Moment4(const Object& obj, const WeightFunction& w) {
    const Grid& grid = obj.grid;
    M[0] = M[1] = M[2] = M[3] = M[4] = 0;
    data_t w_, diff_x, diff_y, val;
    Point<data_t> centroid_wcs = obj.centroid;
    grid.getWCS().transform(centroid_wcs);

    for (long i=0; i< grid.size(); i++) {
      val = obj(i);
      w_ = w(obj.grid(i));
      diff_x = obj.grid(i,0)-centroid_wcs(0);
      diff_y = obj.grid(i,1)-centroid_wcs(1);
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);

      M[0] += pow_int(diff_x,4) * obj(i) * w_;
      M[1] += pow_int(diff_x,3) * diff_y * obj(i) * w_;
      M[2] += pow_int(diff_x,2) * pow_int(diff_y,2) * obj(i) * w_;
      M[3] += diff_x * pow_int(diff_y,3) * obj(i) * w_;
      M[4] += pow_int(diff_y,4) * obj(i) * w_;
    }
  }
  data_t& Moment4::operator()(bool i, bool j, bool k, bool l) {
    return M[i+j+k+l];
  }
  const data_t& Moment4::operator()(bool i, bool j, bool k, bool l) const {
    return M[i+j+k+l];
  }



  Moment6::Moment6() {
    M[0] = M[1] = M[2] = M[3] = M[4] = M[5] = M[6] = 0;
  }
  Moment6::Moment6(const Object& obj, const WeightFunction& w) {
    const Grid& grid = obj.grid;
    M[0] = M[1] = M[2] = M[3] = M[4] = M[5] = M[6] = 0;
    data_t w_, diff_x, diff_y, val;
    Point<data_t> centroid_wcs = obj.centroid;
    grid.getWCS().transform(centroid_wcs);

    for (long i=0; i< grid.size(); i++) {
      val = obj(i);
      w_ = w(obj.grid(i));
      diff_x = obj.grid(i,0)-centroid_wcs(0);
      diff_y = obj.grid(i,1)-centroid_wcs(1);
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);

      M[0] += pow_int(diff_x,6) * obj(i) * w_;
      M[1] += pow_int(diff_x,5) * diff_y * obj(i) * w_;
      M[2] += pow_int(diff_x,4) * pow_int(diff_y,2) * obj(i)* w_;
      M[3] += pow_int(diff_x,3) * pow_int(diff_y,3) * obj(i)* w_;
      M[4] += pow_int(diff_x,2) * pow_int(diff_y,4) * obj(i)* w_;
      M[5] += diff_x * pow_int(diff_y,5) * obj(i) * w_;
      M[6] += pow_int(diff_y,6) * obj(i) * w_;
    }
  }
  data_t& Moment6::operator()(bool i, bool j, bool k, bool l, bool m, bool n) {
    return M[i+j+k+l+m+n];
  }
  const data_t& Moment6::operator()(bool i, bool j, bool k, bool l, bool m, bool n) const {
    return M[i+j+k+l+m+n];
  }

  
  Moment8::Moment8() {
    M[0] = M[1] = M[2] = M[3] = M[4] = M[5] = M[6] = M[7] = M[8] = 0;
  }
  Moment8::Moment8(const Object& obj, const WeightFunction& w) {
    const Grid& grid = obj.grid;
    M[0] = M[1] = M[2] = M[3] = M[4] = M[5] = M[6] = M[7] = M[8] = 0;
    data_t w_, diff_x, diff_y, val;
    Point<data_t> centroid_wcs = obj.centroid;
    grid.getWCS().transform(centroid_wcs);

    for (long i=0; i< grid.size(); i++) {
      val = obj(i);
      w_ = w(obj.grid(i));
      diff_x = obj.grid(i,0)-centroid_wcs(0);
      diff_y = obj.grid(i,1)-centroid_wcs(1);
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);

      M[0] += pow_int(diff_x,8) * val * w_;
      M[1] += pow_int(diff_x,7) * diff_y * val * w_;
      M[2] += pow_int(diff_x,6) * pow_int(diff_y,2) * val * w_;
      M[3] += pow_int(diff_x,5) * pow_int(diff_y,3) * val * w_;
      M[4] += pow_int(diff_x,4) * pow_int(diff_y,4) * val * w_;
      M[5] += pow_int(diff_x,3) * pow_int(diff_y,5) * val * w_;
      M[6] += pow_int(diff_x,2) * pow_int(diff_y,6) * val * w_;
      M[7] += diff_x * pow_int(diff_y,7) * val * w_;
      M[8] += pow_int(diff_y,8) * val * w_;
    }
  }
    data_t& Moment8::operator()(bool i, bool j, bool k, bool l, bool m, bool n, bool o, bool p) {
    return M[i+j+k+l+m+n+o+p];
  }
  const data_t& Moment8::operator()(bool i, bool j, bool k, bool l, bool m, bool n, bool o, bool p) const {
    return M[i+j+k+l+m+n+o+p];
  }




  Moments::Moments() :
    NumVector<data_t>(),
    N(0) { }
  Moments::Moments(int N_) : 
    NumVector<data_t>(pyramid_num(N_+1)),
    N(N_) { }

  Moments& Moments::operator=(const Moments& m) {
    setOrder(m.N);
    NumVector<data_t>::operator=(m);
    return *this;
  }

  Moments::Moments(const Object& obj, const WeightFunction& w, int N_) :
    NumVector<data_t>(pyramid_num(N_+1)),
    N(N_) {
    data_t w_, diff_x, diff_y, val, pow_m, sum_w = 0;
    Point<data_t> centroid_wcs = obj.centroid;
    obj.grid.getWCS().transform(centroid_wcs); // use wcs for centroid
    NumVector<data_t> pow_x(N+1), pow_y(N+1);

    for (long i=0; i< obj.grid.size(); i++) {
      val = obj(i);
      w_ = w(obj.grid(i));
      diff_x = obj.grid(i,0) - centroid_wcs(0);
      diff_y = obj.grid(i,1) - centroid_wcs(1);
      if (obj.weight.size() != 0) {
	w_ *= obj.weight(i);
	sum_w += obj.weight(i);
      }
      for (int i=0; i <= N; i++) {
	pow_x(i) = pow_int(diff_x,i);
	pow_y(i) = pow_int(diff_y,i);
      }
      for(int n=0; n <= N; n++)
	for(int m=0; m <= n; m++) {
	  operator()(m,n-m) += pow_x(m) * pow_y(n-m) * val * w_;
      }
    }
    /*
    if (obj.weight.size() != 0)
      for(int n=0; n <= N; n++)
	for(int m=0; m <= n; m++)
	operator()(m,n-m) /= sum_w;*/
  }

  data_t& Moments::operator()(unsigned int px, unsigned int py) {
    return NumVector<data_t>::operator()(pyramid_num(px+py)+py);
  }

  const data_t& Moments::operator()(unsigned int px, unsigned int py) const {
    return NumVector<data_t>::operator()(pyramid_num(px+py)+py);
  }

  data_t& Moments::operator()(unsigned int i) {
    return NumVector<data_t>::operator()(i);
  }

  const data_t& Moments::operator()(unsigned int i) const {
    return NumVector<data_t>::operator()(i);
  }

  int Moments::getOrder() const {
    return N;
  }

  void Moments::setOrder(int N_) {
    N = N_;
    NumVector<data_t>::resize(pyramid_num(N+1));
  }
  
  int Moments::getIndex(unsigned int px, unsigned int py) const {
    return pyramid_num(px+py)+py;
  }
  unsigned int Moments::pyramid_num(int n) const {
    return (n*(n+1))/2;
  }

} // end namespace

