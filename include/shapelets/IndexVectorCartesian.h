#ifndef INDEXVECTORCARTESIAN_H
#define INDEXVECTORCARTESIAN_H

#include <shapelets/IndexVector.h>
#include <NumMatrix.h>
#include <map>

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
