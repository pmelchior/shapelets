#include "../../include/frame/Moments.h"
#include "../../include/utils/Interpolation.h"

namespace shapelens {

  int SUBPIXEL = 2;
  data_t offset = 1./SUBPIXEL;

  Moment0::Moment0() {
    M = 0;
  }

  Moment0::Moment0(const Object& obj) {
    const Grid& grid = obj.grid;
    M = 0;
    data_t w_;
    for (long i=0; i< grid.size(); i++) {
      //Point<data_t> P = obj.grid(i), P_;
      //for (int n1 = 0; n1 < SUBPIXEL; n1++) {
      // P_(0) = P(0) + (0.5+n1)*offset - 0.5;
      //       for (int n2 = 0; n2 < SUBPIXEL; n2++) {
      // 	P_(1) = P(1) + (0.5+n2)*offset - 0.5;
      // 	w_ = obj.w(P_)/(SUBPIXEL*SUBPIXEL);
      // 	if (obj.weight.size() != 0)
      // 	  w_ *= obj.weight(i);
      // 	M += w_ * obj.interpolate(P_);//Interpolation::bicubic(obj,P_);
      //       }
      //     }
    
      w_ = obj.w(grid(i));
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
  Moment1::Moment1(const Object& obj) {
    const Grid& grid = obj.grid;
    M[0] = M[1] = 0;
    data_t w_;

    for (long i=0; i< grid.size(); i++) {
      // Point<data_t> P = obj.grid(i), P_;
      //     for (int n1 = 0; n1 < SUBPIXEL; n1++) {
      //       P_(0) = P(0) + (0.5+n1)*offset - 0.5;
      //       for (int n2 = 0; n2 < SUBPIXEL; n2++) {
      // 	P_(1) = P(1) + (0.5+n2)*offset - 0.5;
      // 	w_ = obj.w(P_)/(SUBPIXEL*SUBPIXEL);
      // 	if (obj.weight.size() != 0)
      // 	  w_ *= obj.weight(i);
      // 	data_t val = obj.interpolate(P_);//Interpolation::bicubic(obj,P_);
      // 	M[0] += w_ * val * P_(0);
      // 	M[1] += w_ * val * P_(1);
      //       }
      //     }
      w_ = obj.w(grid(i));
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);
      M[0] += obj(i) * grid(i,0) * w_;
      M[1] += obj(i) * grid(i,1) * w_;
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
  Moment2::Moment2(const Object& obj) {
    const Grid& grid = obj.grid;
    M[0] = M[1] = M[2] = 0;
    data_t w_;

    for (long i=0; i< grid.size(); i++) {
      // Point<data_t> P = obj.grid(i), P_;
      //     for (int n1 = 0; n1 < SUBPIXEL; n1++) {
      //       P_(0) = P(0) + (0.5+n1)*offset - 0.5;
      //       for (int n2 = 0; n2 < SUBPIXEL; n2++) {
      // 	P_(1) = P(1) + (0.5+n2)*offset - 0.5;
      // 	w_ = obj.w(P_)/(SUBPIXEL*SUBPIXEL);
      // 	if (obj.weight.size() != 0)
      // 	  w_ *= obj.weight(i);
      // 	data_t val = obj.interpolate(P_);//Interpolation::bicubic(obj,P_);
      // 	M[0] += w_ * val * gsl_pow_2(P_(0)-obj.centroid(0));
      // 	M[1] += w_ * val * (P_(0)-obj.centroid(0))*(P_(1)-obj.centroid(1));
      // 	M[2] += w_ * val * gsl_pow_2(P_(1)-obj.centroid(1));
      //       }
      //     }

      w_ = obj.w(grid(i));
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);
      M[0] += gsl_pow_2(grid(i,0)-obj.centroid(0)) * obj(i) * w_;
      M[1] += (grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
      M[2] += gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
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
  Moment3::Moment3(const Object& obj) {
    const Grid& grid = obj.grid;
    M[0] = M[1] = M[2] = M[3] = 0;
    data_t w_;
  
    for (long i=0; i< grid.size(); i++) {
      w_ = obj.w(grid(i));
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);
      M[0] += gsl_pow_3(grid(i,0)-obj.centroid(0)) * obj(i) * w_;
      M[1] += gsl_pow_2(grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
      M[2] += (grid(i,0)-obj.centroid(0))*gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
      M[3] += gsl_pow_3(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
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
  Moment4::Moment4(const Object& obj) {
    const Grid& grid = obj.grid;
    M[0] = M[1] = M[2] = M[3] = M[4] = 0;
    data_t w_;

    for (long i=0; i< grid.size(); i++) {
      //  Point<data_t> P = obj.grid(i), P_;
      //     for (int n1 = 0; n1 < SUBPIXEL; n1++) {
      //       P_(0) = P(0) + (0.5+n1)*offset - 0.5;
      //       for (int n2 = 0; n2 < SUBPIXEL; n2++) {
      // 	P_(1) = P(1) + (0.5+n2)*offset - 0.5;
      // 	w_ = obj.w(P_)/(SUBPIXEL*SUBPIXEL);
      // 	if (obj.weight.size() != 0)
      // 	  w_ *= obj.weight(i);
      // 	data_t val = obj.interpolate(P_);//Interpolation::bicubic(obj,P_);
      // 	M[0] += w_ * val * gsl_pow_2(P_(0)-obj.centroid(0));
      // 	M[1] += w_ * val * (P_(0)-obj.centroid(0))*(P_(1)-obj.centroid(1));
      // 	M[2] += w_ * val * gsl_pow_2(P_(1)-obj.centroid(1));
      //       }
      //     }
    
      w_ = obj.w(grid(i));
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);
      M[0] += gsl_pow_4(grid(i,0)-obj.centroid(0)) * obj(i) * w_;
      M[1] += gsl_pow_3(grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
      M[2] += gsl_pow_2(grid(i,0)-obj.centroid(0))*gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
      M[3] += (grid(i,0)-obj.centroid(0))*gsl_pow_3(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
      M[4] += gsl_pow_4(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
    }
  }
  data_t& Moment4::operator()(bool i, bool j, bool k, bool l) {
    return M[i+j+k+l];
  }
  const data_t& Moment4::operator()(bool i, bool j, bool k, bool l) const {
    return M[i+j+k+l];
  }



  Moment5::Moment5() {
    M[0] = M[1] = M[2] = M[3] = M[4] = M[5] = 0;
  }
  Moment5::Moment5(const Object& obj) {
    const Grid& grid = obj.grid;
    M[0] = M[1] = M[2] = M[3] = M[4] = M[5] = 0;
    data_t w_;
  
    for (long i=0; i< grid.size(); i++) {
      w_ = obj.w(grid(i));
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);
      M[0] += gsl_pow_5(grid(i,0)-obj.centroid(0)) * obj(i) * w_;
      M[1] += gsl_pow_4(grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
      M[2] += gsl_pow_3(grid(i,0)-obj.centroid(0))*gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i)* w_;
      M[3] += gsl_pow_2(grid(i,0)-obj.centroid(0))*gsl_pow_3(grid(i,1)-obj.centroid(1)) * obj(i)* w_;
      M[4] += (grid(i,0)-obj.centroid(0))*gsl_pow_4(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
      M[5] += gsl_pow_5(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
    }
  }
  data_t& Moment5::operator()(bool i, bool j, bool k, bool l, bool m) {
    return M[i+j+k+l+m];
  }
  const data_t& Moment5::operator()(bool i, bool j, bool k, bool l, bool m) const {
    return M[i+j+k+l+m];
  }



  Moment6::Moment6() {
    M[0] = M[1] = M[2] = M[3] = M[4] = M[5] = M[6] = 0;
  }
  Moment6::Moment6(const Object& obj) {
    const Grid& grid = obj.grid;
    M[0] = M[1] = M[2] = M[3] = M[4] = M[5] = M[6] = 0;
    data_t w_;

    for (long i=0; i< grid.size(); i++) {
      w_ = obj.w(grid(i));
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);
      M[0] += gsl_pow_6(grid(i,0)-obj.centroid(0)) * obj(i) * w_;
      M[1] += gsl_pow_5(grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
      M[2] += gsl_pow_4(grid(i,0)-obj.centroid(0))*gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i)* w_;
      M[3] += gsl_pow_3(grid(i,0)-obj.centroid(0))*gsl_pow_3(grid(i,1)-obj.centroid(1)) * obj(i)* w_;
      M[4] += gsl_pow_2(grid(i,0)-obj.centroid(0))*gsl_pow_4(grid(i,1)-obj.centroid(1)) * obj(i)* w_;
      M[5] += (grid(i,0)-obj.centroid(0))*gsl_pow_5(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
      M[6] += gsl_pow_6(grid(i,1)-obj.centroid(1)) * obj(i) * w_;
    }
  }
  data_t& Moment6::operator()(bool i, bool j, bool k, bool l, bool m, bool n) {
    return M[i+j+k+l+m+n];
  }
  const data_t& Moment6::operator()(bool i, bool j, bool k, bool l, bool m, bool n) const {
    return M[i+j+k+l+m+n];
  }




  MomentsOrdered::MomentsOrdered() :
    NumVector<data_t>(),
    N(0) { }
  MomentsOrdered::MomentsOrdered(int N_) : 
    NumVector<data_t>(pyramid_num(N_+1)),
    N(N_) { }
  MomentsOrdered::MomentsOrdered(const Object& obj, int N_) :
    NumVector<data_t>(pyramid_num(N_+1)),
    N(N_) {
    data_t w_, diff_x, diff_y, val;
    for (long i=0; i< obj.grid.size(); i++) {
      val = obj(i);
      w_ = obj.w(obj.grid(i));
      diff_x = obj.grid(i,0)-obj.centroid(0);
      diff_y = obj.grid(i,1)-obj.centroid(1);
      if (obj.weight.size() != 0)
	w_ *= obj.weight(i);
      for(int n=0; n <= N; n++)
	for(int m=0; m <= n; m++)
	  operator()(m,n-m) += gsl_pow_int(diff_x,m) * gsl_pow_int(diff_y,n-m) * val * w_;
    }
  }
  data_t& MomentsOrdered::operator()(unsigned int px, unsigned int py) {
    return NumVector<data_t>::operator()(pyramid_num(px+py)+py);
  }
  const data_t& MomentsOrdered::operator()(unsigned int px, unsigned int py) const {
    return NumVector<data_t>::operator()(pyramid_num(px+py)+py);
  }
  int MomentsOrdered::getOrder() const {
    return N;
  }
  int MomentsOrdered::getIndex(unsigned int px, unsigned int py) const {
    return pyramid_num(px+py)+py;
  }
  unsigned int MomentsOrdered::pyramid_num(int n) const {
    return (n*(n+1))/2;
  }

} // end namespace

