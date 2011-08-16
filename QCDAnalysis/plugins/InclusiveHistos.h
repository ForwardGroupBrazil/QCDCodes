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
#include "TProfile.h"

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
    std::string mFileName,mTreeName,mDirName,mCaloJECres,mPFJECres,mCaloYBiasCor,mPFYBiasCor,mPUHistName,mPUFileName;
    std::vector<double> mYBND,mPTBND;
    std::vector<std::string> mTriggers;
    std::vector<int> mTrigIndex;    

    edm::Service<TFileService> fs;
    TTree *mTree; 
    TFile *mInf,*mPUf;
    TDirectoryFile *mDir;
    TH1F *mBSx,*mBSy,*mBSz,*mNPV,*mPVx,*mPVy,*mPVz,*mPUh;
    TH1F *mNPFJets[100][6],*mNCaloJets[100][6];
    TH1F *mPFJetMulti[100],*mCaloJetMulti[100],*mPFMETovSUMET[100][6],*mCaloMETovSUMET[100][6];
    TH1F *mGenPt[6],*mGenX[6],*mPFPt[100][6],*mPFNormPt[100][6],*mPFX[100][6],*mPFNormX[100][6],
         *mCaloPt[100][6],*mCaloNormPt[100][6],*mCaloX[100][6],*mCaloNormX[100][6],
         *mCHF[100][6],*mNHF[100][6],*mPHF[100][6],*mELF[100][6],
         *mN90hits[100][6],*mEMF[100][6],*mNTrkCalo[100][6],*mNTrkVtx[100][6],*mfHPD[100][6];
    TProfile *mEMFVsRun[100][6],*mNTrkCaloVsRun[100][6],*mNTrkVtxVsRun[100][6];
    TProfile *mNHFVsRun[100][6],*mPHFVsRun[100][6],*mCHFVsRun[100][6],*mELFVsRun[100][6];
    TProfile *mCaloRhoVsRun,*mPFRhoVsRun,*mCaloRhoVsNPV,*mPFRhoVsNPV; 
    //---- TREE variable --------
    QCDEvent *mEvent;
    
};

#endif
