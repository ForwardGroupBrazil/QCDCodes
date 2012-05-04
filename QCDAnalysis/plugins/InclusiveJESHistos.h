#ifndef InclusiveJESHistos_h
#define InclusiveJESHistos_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "KKousour/QCDAnalysis/interface/QCDJet.h"
#include "KKousour/QCDAnalysis/interface/QCDEvent.h"
#include "KKousour/QCDAnalysis/interface/QCDEventHdr.h"
#include "KKousour/QCDAnalysis/interface/QCDCaloJet.h"
#include "KKousour/QCDAnalysis/interface/QCDPFJet.h"
#include "KKousour/QCDAnalysis/interface/QCDMET.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"

class InclusiveJESHistos : public edm::EDAnalyzer 
{
  public:
    explicit InclusiveJESHistos(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~InclusiveJESHistos();

  private:  
    int getBin(double x, const std::vector<double>& boundaries); 
    int    mNEvents;
    double mMaxMETovSumET;
    std::vector<double> mMinPt;
    std::string mTreeName,mDirName;
    std::vector<double> mYBND,mPTBND;
    std::vector<std::string> mTriggers;
    std::vector<std::string> mFileNames;
    std::vector<std::string> mJECUncSrcNames;
    std::vector<std::vector<int> > mTrigIndex;    

    edm::Service<TFileService> fs;
    TH1F *mPt[100][6],*mPtUp[100][6][50],*mPtDown[100][6][50];
    TH1F *mPtTotal[6],*mPtUpTotal[6][50],*mPtDownTotal[6][50];
    TH1F *mPtRatioUpTotal[6][50],*mPtRatioDownTotal[6][50];
    //---- TREE variable --------
    QCDEvent *mEvent;
    
};

#endif
