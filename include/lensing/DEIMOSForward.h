#ifndef SHAPELENS_DEIMOSFORWARD_H
#define SHAPELENS_DEIMOSFORWARD_H
#include "DEIMOS.h"
#include <vector>

namespace shapelens {
  typedef std::vector<Object> MultiExposureObject;
  typedef std::vector<Moments> MultiExposureMoments;
  class DEIMOSForward : public DEIMOS {
  public:
    DEIMOSForward(const MultiExposureObject& meo, const MultiExposureMoments& mepsf, int N, int C, data_t flux, data_t width);
    data_t operator()(const NumVector<data_t>& p);
  protected:
    void convolveExposure(unsigned int k);
    int K;
  private:
    void computeMomentsFromGuess();
    const MultiExposureObject& meo;
    const MultiExposureMoments& mepsf;
    MultiExposureMoments mem_, mem;
    Moments mo0;
    data_t width;
    std::vector<NumMatrix<data_t> > meP, meS;
  };
} // end namespace
#endif
