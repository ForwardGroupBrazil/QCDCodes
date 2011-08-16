#ifndef ResponseHistos_h
#define ResponseHistos_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "KKousour/QCDAnalysis/interface/QCDEvent.h"
#include "TTree.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"

class ResponseHistos : public edm::EDAnalyzer 
{
  public:
    explicit ResponseHistos(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~ResponseHistos();

  private:  
    int getBin(double x, const std::vector<double>& boundaries);
    //---- configurable parameters --------   
    std::string mFileName,mTreeName,mDirName,mPFBiasCorName;
    std::vector<double> mPtBND;
    std::vector<double> mYBND;
    std::vector<double> mFineYBND;
    double mMaxDR;
    int    mNEvents;
    int    mMaxJets;

    edm::Service<TFileService> fs;
    QCDEvent *mEvent;
    TTree *mTree; 
    TFile *mInf;
    TDirectoryFile *mDir;
    TH2F *mCaloRspVsY[20],*mCaloRspVsPt[20],*mPFRspVsY[20],*mPFRspVsPt[20],*mPFDYVsPt[20],*mCaloDYVsPt[20];
    TProfile *mCaloDYVsY[20],*mPFDYVsY[20],*mCaloDYVsAbsY[20],*mPFDYVsAbsY[20],*mCaloPtVsAbsY[20],*mPFPtVsAbsY[20];
    TProfile *mPFDCorYVsY[20];
    TH2F *mGenYVsPt;
};

#endif
