#ifndef SHAPELENS_DEIMOSFORWARD_H
#define SHAPELENS_DEIMOSFORWARD_H
#include "DEIMOS.h"
#include <vector>

namespace shapelens {
  typedef std::vector<Object> MultiExposureObject;
  typedef std::vector<Moments> MultiExposureMoments;
  class DEIMOSForward : public DEIMOS {
  public:
    DEIMOSForward(const MultiExposureObject& meo, const MultiExposureObject& mepsf, int N, int C, data_t flux, data_t width);
    data_t operator()(const NumVector<data_t>& p);
    std::vector<DEIMOS> meD;
  private:
    void computeMomentsFromGuess();
    data_t getWeightFunctionScale(const Moments& m) const;
    const MultiExposureObject& meo;
    const MultiExposureObject& mepsf;
    MultiExposureMoments mem;
    Moments mo0;
    std::vector<NumMatrix<data_t> > meP;
    std::vector<DEIMOS::PSFMultiScale> mePSFMultiScale;
  protected:
    void convolveExposure(unsigned int k);
    int K;
  };
} // end namespace
#endif
