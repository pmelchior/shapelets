#include "../../include/frame/Moments.h"

using namespace shapelens;

Moment0::Moment0() {
  M = 0;
}

Moment0::Moment0(const Object& obj) {
  const Grid& grid = obj.grid;
  M = 0;
  // check if weights are available: if yes, use them
  if (obj.weight.size() != 0) {
    data_t sum_weights = 0;
    for (long i=0; i< grid.size(); i++) {
      M += obj(i) * obj.weight(i);
      sum_weights += obj.weight(i);
    }
    if (sum_weights != 0) 
      M /= sum_weights;
  }
  else // unweigthed
    for (long i=0; i< grid.size(); i++)
      M += obj(i);
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
  // check if weights are available: if yes, use them
  if (obj.weight.size() != 0) {
    data_t sum_weights = 0;
    for (long i=0; i< grid.size(); i++) {
      M[0] += obj(i) * grid(i,0) * obj.weight(i);
      M[1] += obj(i) * grid(i,1) * obj.weight(i);
      sum_weights += obj.weight(i);
    }
    if (sum_weights != 0) {
      M[0] /= sum_weights;
      M[1] /= sum_weights;
    }
  }
  else { // unweighted
    for (long i=0; i< grid.size(); i++) {
      M[0] += grid(i,0) * obj(i);
      M[1] += grid(i,1) * obj(i);
    }
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
  // check if weights are available: if yes, use them
  if (obj.weight.size() != 0) {
    data_t sum_weights = 0;
    for (long i=0; i< grid.size(); i++) {
      M[0] += gsl_pow_2(grid(i,0)-obj.centroid(0)) * obj(i) * obj.weight(i);
      M[1] += (grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[2] += gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      sum_weights += obj.weight(i);
    }
    if (sum_weights != 0) {
      M[0] /= sum_weights;
      M[1] /= sum_weights;
      M[2] /= sum_weights;
    }
  }
  else { // unweighted
    for (long i=0; i< grid.size(); i++) {
      M[0] += gsl_pow_2(grid(i,0)-obj.centroid(0)) * obj(i);
      M[1] += (grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i);
      M[2] += gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i);
    }
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
  // check if weights are available: if yes, use them
  if (obj.weight.size() != 0) {
    data_t sum_weights = 0;
    for (long i=0; i< grid.size(); i++) {
      M[0] += gsl_pow_3(grid(i,0)-obj.centroid(0)) * obj(i) * obj.weight(i);
      M[1] += gsl_pow_2(grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[2] += (grid(i,0)-obj.centroid(0))*gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[3] += gsl_pow_3(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      sum_weights += obj.weight(i);
    }
    if (sum_weights != 0) {
      M[0] /= sum_weights;
      M[1] /= sum_weights;
      M[2] /= sum_weights;
      M[3] /= sum_weights;
    }
  }
  else { // unweighted
    for (long i=0; i< grid.size(); i++) {
      M[0] += gsl_pow_3(grid(i,0)-obj.centroid(0)) * obj(i);
      M[1] += gsl_pow_2(grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i);
      M[2] += (grid(i,0)-obj.centroid(0))*gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i);
      M[3] += gsl_pow_3(grid(i,1)-obj.centroid(1)) * obj(i);
    }
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
  // check if weights are available: if yes, use them
  if (obj.weight.size() != 0) {
    data_t sum_weights = 0;
    for (long i=0; i< grid.size(); i++) {
      M[0] += gsl_pow_4(grid(i,0)-obj.centroid(0)) * obj(i) * obj.weight(i);
      M[1] += gsl_pow_3(grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[2] += gsl_pow_2(grid(i,0)-obj.centroid(0))*gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[3] += (grid(i,0)-obj.centroid(0))*gsl_pow_3(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[4] += gsl_pow_4(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      sum_weights += obj.weight(i);
    }
    if (sum_weights != 0) {
      M[0] /= sum_weights;
      M[1] /= sum_weights;
      M[2] /= sum_weights;
      M[3] /= sum_weights;
      M[4] /= sum_weights;
    }
  }
  else { // unweighted
    for (long i=0; i< grid.size(); i++) {
      M[0] += gsl_pow_4(grid(i,0)-obj.centroid(0)) * obj(i);
      M[1] += gsl_pow_3(grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i);
      M[2] += gsl_pow_2(grid(i,0)-obj.centroid(0))*gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i);
      M[3] += (grid(i,0)-obj.centroid(0))*gsl_pow_3(grid(i,1)-obj.centroid(1)) * obj(i);
      M[4] += gsl_pow_4(grid(i,1)-obj.centroid(1)) * obj(i);
    }
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
  // check if weights are available: if yes, use them
  if (obj.weight.size() != 0) {
    data_t sum_weights = 0;
    for (long i=0; i< grid.size(); i++) {
      M[0] += gsl_pow_5(grid(i,0)-obj.centroid(0)) * obj(i) * obj.weight(i);
      M[1] += gsl_pow_4(grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[2] += gsl_pow_3(grid(i,0)-obj.centroid(0))*gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[3] += gsl_pow_2(grid(i,0)-obj.centroid(0))*gsl_pow_3(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[4] += (grid(i,0)-obj.centroid(0))*gsl_pow_4(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[5] += gsl_pow_5(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      sum_weights += obj.weight(i);
    }
    if (sum_weights != 0) {
      M[0] /= sum_weights;
      M[1] /= sum_weights;
      M[2] /= sum_weights;
      M[3] /= sum_weights;
      M[4] /= sum_weights;
      M[5] /= sum_weights;
    }
  }
  else { // unweighted
    for (long i=0; i< grid.size(); i++) {
      M[0] += gsl_pow_5(grid(i,0)-obj.centroid(0)) * obj(i);
      M[1] += gsl_pow_4(grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i);
      M[2] += gsl_pow_3(grid(i,0)-obj.centroid(0))*gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i);
      M[3] += gsl_pow_2(grid(i,0)-obj.centroid(0))*gsl_pow_3(grid(i,1)-obj.centroid(1)) * obj(i);
      M[4] += (grid(i,0)-obj.centroid(0))*gsl_pow_4(grid(i,1)-obj.centroid(1)) * obj(i);
      M[5] += gsl_pow_5(grid(i,1)-obj.centroid(1)) * obj(i);
    }
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
  // check if weights are available: if yes, use them
  if (obj.weight.size() != 0) {
    data_t sum_weights = 0;
    for (long i=0; i< grid.size(); i++) {
      M[0] += gsl_pow_6(grid(i,0)-obj.centroid(0)) * obj(i) * obj.weight(i);
      M[1] += gsl_pow_5(grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[2] += gsl_pow_4(grid(i,0)-obj.centroid(0))*gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[3] += gsl_pow_3(grid(i,0)-obj.centroid(0))*gsl_pow_3(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[4] += gsl_pow_2(grid(i,0)-obj.centroid(0))*gsl_pow_4(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[5] += (grid(i,0)-obj.centroid(0))*gsl_pow_5(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      M[6] += gsl_pow_6(grid(i,1)-obj.centroid(1)) * obj(i) * obj.weight(i);
      sum_weights += obj.weight(i);
    }
    if (sum_weights != 0) {
      M[0] /= sum_weights;
      M[1] /= sum_weights;
      M[2] /= sum_weights;
      M[3] /= sum_weights;
      M[4] /= sum_weights;
      M[5] /= sum_weights;
      M[6] /= sum_weights;
    }
  }
  else { // unweighted
    for (long i=0; i< grid.size(); i++) {
      M[0] += gsl_pow_6(grid(i,0)-obj.centroid(0)) * obj(i);
      M[1] += gsl_pow_5(grid(i,0)-obj.centroid(0))*(grid(i,1)-obj.centroid(1)) * obj(i);
      M[2] += gsl_pow_4(grid(i,0)-obj.centroid(0))*gsl_pow_2(grid(i,1)-obj.centroid(1)) * obj(i);
      M[3] += gsl_pow_3(grid(i,0)-obj.centroid(0))*gsl_pow_3(grid(i,1)-obj.centroid(1)) * obj(i);
      M[4] += gsl_pow_2(grid(i,0)-obj.centroid(0))*gsl_pow_4(grid(i,1)-obj.centroid(1)) * obj(i);
      M[5] += (grid(i,0)-obj.centroid(0))*gsl_pow_5(grid(i,1)-obj.centroid(1)) * obj(i);
      M[6] += gsl_pow_6(grid(i,1)-obj.centroid(1)) * obj(i);
    }
  }
}
data_t& Moment6::operator()(bool i, bool j, bool k, bool l, bool m, bool n) {
  return M[i+j+k+l+m+n];
}
const data_t& Moment6::operator()(bool i, bool j, bool k, bool l, bool m, bool n) const {
  return M[i+j+k+l+m+n];
}

