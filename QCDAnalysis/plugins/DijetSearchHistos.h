#ifndef DijetSearchHistos_h
#define DijetSearchHistos_h

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

class DijetSearchHistos : public edm::EDAnalyzer 
{
  public:
    explicit DijetSearchHistos(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~DijetSearchHistos();

  private:  
    int findRun(int x, const std::vector<int>& runs);
    //---- configurable parameters --------   
    bool   mIsMC;
    int    mJetID;
    int    mHCALNoise;
    int    mNEvents;
    double mMinPt;
    double mMinMass;
    double mMaxEta;
    double mMaxDeta;
    double mChiIN;
    double mChiOUT; 
    std::string mFileName,mTreeName,mDirName;
    std::vector<double> mMassBND;
    std::vector<std::string> mTriggers;
    std::vector<int> mTrigIndex;    

    edm::Service<TFileService> fs;
    TTree *mTree; 
    TFile *mInf;
    TDirectoryFile *mDir;
    TH1F *mPFMETovSUMET,*mCaloMETovSUMET;
    TH1F *mPFPt,*mCaloPt,*mPFMass,*mCaloMass,*mPFEta,*mCaloEta,*mPFDeta,*mCaloDeta,*mPFDphi,*mCaloDphi,*mPFPtJJ,*mCaloPtJJ,*mPFEtaBoost,*mCaloEtaBoost,
         *mCHF,*mNHF,*mPHF,*mPFMassIN,*mPFMassOUT,*mCaloMassIN,*mCaloMassOUT,
         *mPFChi,*mCaloChi,*mPFChiIN,*mCaloChiIN,*mPFChiOUT,*mCaloChiOUT,
         *mN90hits,*mEMF,*mNTrkCalo,*mNTrkVtx,*mfHPD;
    TProfile *mEMFVsRun,*mNTrkCaloVsRun,*mNTrkVtxVsRun;
    TProfile *mNHFVsRun,*mPHFVsRun,*mCHFVsRun; 
    //---- TREE variable --------
    QCDEvent *mEvent;
    
};

#endif
