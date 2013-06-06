#ifndef QGLCalculator_h
#define QGLCalculator_h

#include "TMVA/Reader.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"

class QGLCalculator
{
  public:
    QGLCalculator();
    ~QGLCalculator();
    float getQGL(cmg::PFJet const& jet,float rho);
    
  private:
    int FindIndex(int N,const float BND[],float x);
    std::string PT_CAT_[7];
    std::string ETA_CAT_[3];
    float PT_BND_[8],ETA_BND_[4];
    float COR_RHO_[7][3][6];
    TMVA::Reader *readerQGL_[7][3];
    float varQGL_[5];
};
#endif
