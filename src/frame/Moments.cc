#include "../../include/frame/Moments.h"

using namespace shapelens;

Quadrupole::Quadrupole() {
  m[0] = m[1] = m[2] = 0;
}
data_t& Quadrupole::operator()(bool i, bool j) {
  // count number of ones in index set
  int n = i+j;
  return m[n];
}
const data_t& Quadrupole::operator()(bool i, bool j) const {
// count number of ones in index set
  int n = i+j;
  return m[n];
}

Octupole::Octupole() {
  m[0] = m[1] = m[2] = m[3] = 0;
}
data_t& Octupole::operator()(bool i, bool j, bool k) {
  // count number of ones in index set
  int n = i+j+k;
  return m[n];
}
const data_t& Octupole::operator()(bool i, bool j, bool k) const {
  // count number of ones in index set
  int n = i+j+k;
  return m[n];
}

Hexadecupole::Hexadecupole() {
  m[0] = m[1] = m[2] = m[3] = m[4] = 0;
}
data_t& Hexadecupole::operator()(bool i, bool j, bool k, bool l) {
  // count number of ones in index set
  int n = i+j+k+l;
  return m[n];
}
const data_t& Hexadecupole::operator()(bool i, bool j, bool k, bool l) const {
  // count number of ones in index set
  int n = i+j+k+l;
  return m[n];
}
