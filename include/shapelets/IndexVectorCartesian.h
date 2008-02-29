#ifndef INDEXVECTORCARTESIAN_H
#define INDEXVECTORCARTESIAN_H

#include <shapelets/IndexVector.h>
#include <NumMatrix.h>
#include <map>

class IndexVectorCartesian : public IndexVector {
 public:
  IndexVectorCartesian();
  IndexVectorCartesian(unsigned int nmax);
  virtual ~IndexVectorCartesian() {};
  virtual unsigned int getNMax() const;
  virtual unsigned int getNCoeffs() const;
  virtual void setNMax(unsigned int nmax);
  virtual int getState1(unsigned int index) const;
  virtual int getState2(unsigned int index) const;
  virtual unsigned int getIndex1(unsigned int index) const;
  virtual unsigned int getIndex2(unsigned int index) const;
  unsigned int getIndex(int n1, int n2) const;
 private:
  unsigned int nmax;
  std::map<unsigned int, std::pair<unsigned int, unsigned int> > cartesian;
  NumMatrix<unsigned int> cMatrix;
  void computeIndexMaps();
};

#endif
