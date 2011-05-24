#ifndef InclusiveHistos_h
#define InclusiveHistos_h

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
#include "TH1F.h"
#include "TFile.h"

class InclusiveHistos : public edm::EDAnalyzer 
{
  public:
    explicit InclusiveHistos(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~InclusiveHistos();

  private:  
    int getBin(double x, const std::vector<double>& boundaries); 
    int findRun(int x, const std::vector<int>& runs);
    //---- configurable parameters --------   
    bool   mIsMC;
    int    mJetID;
    int    mHCALNoise;
    int    mNEvents;
    std::vector<double> mMinPt;
    std::string mFileName,mTreeName,mDirName;
    std::vector<double> mYBND,mPTBND;
    std::vector<std::string> mTriggers;
    std::vector<int> mTrigIndex;    

    edm::Service<TFileService> fs;
    TTree *mTree; 
    TFile *mInf;
    TDirectoryFile *mDir;
    TH1F *mNPFJets[30][10],*mNCaloJets[30][10];
    TH1F *mPFJetMulti[30],*mCaloJetMulti[30],*mPFMETovSUMET[30][10],*mCaloMETovSUMET[30][10];
    TH1F *mPFPt[30][10],*mPFNormPt[30][10],*mPFX[30][10],*mPFNormX[30][10],
         *mCaloPt[30][10],*mCaloNormPt[30][10],*mCaloX[30][10],*mCaloNormX[30][10],
         *mCHF[30][10],*mNHF[30][10],*mPHF[30][10],
         *mN90hits[30][10],*mEMF[30][10],*mNTrkCalo[30][10],*mNTrkVtx[30][10],*mfHPD[30][10];
    //---- TREE variable --------
    QCDEvent *mEvent;
    
};

#endif
