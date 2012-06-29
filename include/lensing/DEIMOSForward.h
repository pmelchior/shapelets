#ifndef SHAPELENS_DEIMOSFORWARD_H
#define SHAPELENS_DEIMOSFORWARD_H
#include "DEIMOSElliptical.h"
#include <vector>

namespace shapelens {
  typedef std::vector<Object> MultiExposureObject;
  typedef std::vector<Moments> MultiExposureMoments;
  class DEIMOSForward : public DEIMOSElliptical {
  public:
    DEIMOSForward(const MultiExposureObject& meo, const MultiExposureObject& mepsf, int N, int C, data_t fiducial_width);
    DEIMOSForward(const MultiExposureObject& meo, const std::vector<DEIMOS::PSFMultiScale>& mePSFMultiScale, int N, int C, data_t fiducial_width);
    std::vector<DEIMOSElliptical> meD;
  private:
    void computeMomentsFromGuess();
    data_t getWeightFunctionScale(unsigned int k) const;
    const MultiExposureObject& meo;
    const MultiExposureObject& mepsf;
    MultiExposureMoments mem;
    Moments mo0;
    std::vector<NumMatrix<data_t> > meP;
    std::vector<DEIMOS::PSFMultiScale> mePSFMultiScale;
    data_t fiducial_width;
    std::vector<data_t> meSaveScale, meTroubleScale;
  protected:
    void initialize();
    void minimize();
    void convolveExposure(unsigned int k);
    int K;
  };
} // end namespace
#endif
