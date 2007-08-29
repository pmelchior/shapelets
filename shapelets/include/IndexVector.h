#ifndef INDEXVECTOR_H
#define INDEXVECTOR_H

#include <NumMatrix.h>

class IndexVector : public NumMatrix<int> {
 public:
  IndexVector();
  IndexVector (int nmax);
  void setOrder(int nmax); 
  int getOrder() const;
  int getNCoeffs() const;
  int getN1(unsigned int i) const;
  int getN2(unsigned int i) const;
  int getN(unsigned int i) const;
 private:
  void computeIndexVector();
  unsigned int nmax;
  NumMatrix<int> nVector;
};
#endif
