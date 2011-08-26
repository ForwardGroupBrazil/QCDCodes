#ifndef MultijetSearchHistos_h
#define MultijetSearchHistos_h

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
#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"

class MultijetSearchHistos : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit MultijetSearchHistos(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~MultijetSearchHistos();

  private:  
    int getBin(double x, const std::vector<double>& boundaries);
    //---- configurable parameters --------   
    bool   mIsMC;
    int    mJetID;
    int    mHCALNoise;
    int    mNEvents;
    int    mRank;
    double mScale4J;
    double mOffset4J;
    double mMaxEta;
    double mMinHT;
    std::string mFileName,mTreeName,mDirName;
    std::vector<double> mMinPt;
    std::vector<double> mPtHatLumi,mPtHatBnd;
    std::vector<std::string> mTriggers;
    std::vector<int> mTrigIndex;    

    edm::Service<TFileService> fs;
    TTree *mTree; 
    TFile *mInf;
    TDirectoryFile *mDir;
    TH1F *mMETovSUMET,*mHT,*mHTAll,*mHT4J,*mHT4Jcut,*mMAllJ,*mM4J,*mM4Jcut,*mM2J,*mDR,*mJetMulti,*mPtHat,*mPtHatAll,*mPtRatio;
    TH2F *mM4JvsHT4J;
    TH1F *mPt[8],*mPhi[8],*mEta[8],*mCHF[8],*mNHF[8],*mPHF[8],*mELF[8],*mMUF[8],*mBeta[8];
    //---- TREE variable --------
    QCDEvent *mEvent;
};

#endif
