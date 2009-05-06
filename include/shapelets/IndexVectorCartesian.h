#ifndef SHAPELENS_INDEXVECTORCARTESIAN_H
#define SHAPELENS_INDEXVECTORCARTESIAN_H

#include <map>
#include <numla/NumMatrix.h>
#include "IndexVector.h"

namespace shapelens {
class IndexVectorCartesian : public IndexVector {
 public:
  IndexVectorCartesian();
  IndexVectorCartesian(unsigned int nmax);
  virtual ~IndexVectorCartesian() {};
  virtual int getNMax() const;
  virtual int getNCoeffs() const;
  virtual void setNMax(unsigned int nmax);
  virtual int getState1(unsigned int index) const;
  virtual int getState2(unsigned int index) const;
  virtual int getIndex1(unsigned int index) const;
  virtual int getIndex2(unsigned int index) const;
  virtual int getIndex(int n1, int n2) const;
 private:
  unsigned int nmax;
  std::map<unsigned int, std::pair<unsigned int, unsigned int> > cartesian;
  NumMatrix<unsigned int> cMatrix;
  void computeIndexMaps();
};
} // end namespace
#endif
