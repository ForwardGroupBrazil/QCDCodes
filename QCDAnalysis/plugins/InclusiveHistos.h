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
    TH1F *mNPFJets[50][6],*mNCaloJets[50][6];
    TH1F *mPFJetMulti[50],*mCaloJetMulti[50],*mPFMETovSUMET[50][6],*mCaloMETovSUMET[50][6];
    TH1F *mPFPt[50][6],*mPFNormPt[50][6],*mPFX[50][6],*mPFNormX[50][6],
         *mCaloPt[50][6],*mCaloNormPt[50][6],*mCaloX[50][6],*mCaloNormX[50][6],
         *mCHF[50][6],*mNHF[50][6],*mPHF[50][6],
         *mN90hits[50][6],*mEMF[50][6],*mNTrkCalo[50][6],*mNTrkVtx[50][6],*mfHPD[50][6];
    //---- TREE variable --------
    QCDEvent *mEvent;
    
};

#endif
