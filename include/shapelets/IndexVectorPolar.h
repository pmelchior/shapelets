#ifndef INDEXVECTORPOLAR_H
#define INDEXVECTORPOLAR_H

#include <shapelets/IndexVector.h>
#include <NumMatrix.h>
#include <map>

class IndexVectorPolar : public IndexVector {
 public:
  IndexVectorPolar();
  IndexVectorPolar(unsigned int nmax);
  virtual ~IndexVectorPolar() {};
  virtual unsigned int getNMax() const;
  virtual unsigned int getNCoeffs() const;
  virtual void setNMax(unsigned int nmax);
  virtual int getState1(unsigned int index) const;
  virtual int getState2(unsigned int index) const;
  virtual unsigned int getIndex1(unsigned int index) const;
  virtual unsigned int getIndex2(unsigned int index) const;
  virtual unsigned int getIndex(int n, int m) const;
 private:
  unsigned int nmax;
  std::map<unsigned int, std::pair<unsigned int, int> > polar;
  NumMatrix<unsigned int> pMatrix;
  int mIndex(unsigned int n, int m) const;
  void computeIndexMaps();
};

#endif
