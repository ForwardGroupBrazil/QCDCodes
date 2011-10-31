#ifndef MultijetSearchGENHistos_h
#define MultijetSearchGENHistos_h

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

class MultijetSearchGENHistos : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit MultijetSearchGENHistos(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~MultijetSearchGENHistos();

  private:  
    int getBin(double x, const std::vector<double>& boundaries);
    //---- configurable parameters --------   
    int    mNEvents;
    int    mRank;
    double mMaxEta;
    double mMinHT;
    double mScale4J;
    double mOffset4J;
    std::string mFileName,mTreeName,mDirName;
    std::vector<double> mMinPt;
    std::vector<double> mPtHatLumi,mPtHatBnd;

    edm::Service<TFileService> fs;
    TTree *mTree; 
    TFile *mInf;
    TDirectoryFile *mDir;
    TH1F *mHT,*mHTAll,*mHT4J,*mHT4Jcut,*mMAllJ,*mM4J,*mM4Jcut,*mM2J,*mDR,*mJetMulti,*mPtHat,*mPtHatAll,*mPtRatio;
    TH2F *mM4JvsHT4J;
    TH1F *mPt[8],*mPhi[8],*mEta[8];
    //---- TREE variable --------
    QCDEvent *mEvent;
};

#endif
