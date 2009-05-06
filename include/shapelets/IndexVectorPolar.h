#ifndef SHAPELENS_INDEXVECTORPOLAR_H
#define SHAPELENS_INDEXVECTORPOLAR_H

#include <map>
#include <numla/NumMatrix.h>
#include "IndexVector.h"

namespace shapelens {
class IndexVectorPolar : public IndexVector {
 public:
  IndexVectorPolar();
  IndexVectorPolar(unsigned int nmax);
  virtual ~IndexVectorPolar() {};
  virtual int getNMax() const;
  virtual int getNCoeffs() const;
  virtual void setNMax(unsigned int nmax);
  virtual int getState1(unsigned int index) const;
  virtual int getState2(unsigned int index) const;
  virtual int getIndex1(unsigned int index) const;
  virtual int getIndex2(unsigned int index) const;
  virtual int getIndex(int n, int m) const;
 private:
  unsigned int nmax;
  std::map<unsigned int, std::pair<unsigned int, int> > polar;
  NumMatrix<unsigned int> pMatrix;
  int mIndex(unsigned int n, int m) const;
  void computeIndexMaps();
};
} // end namespace
#endif
