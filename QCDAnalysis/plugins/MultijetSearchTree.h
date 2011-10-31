#ifndef MultijetSearchTree_h
#define MultijetSearchTree_h

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
#include "TFile.h"

class MultijetSearchTree : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit MultijetSearchTree(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~MultijetSearchTree();

  private:  
    std::vector<int> findIndices(const std::vector<int>& v, int rank);
    int getBin(double x, const std::vector<double>& boundaries);
    //---- configurable parameters --------   
    bool   mIsMC;
    bool   mIsPreScaled;
    int    mJetID;
    int    mHCALNoise;
    int    mNEvents;
    double mEtaMax;
    double mPtMin;
    double mPUTagMin;
    std::string mTreeName,mDirName;
    std::vector<double> mPtHatLumi,mPtHatBnd;
    std::vector<std::string> mTriggers,mFileNames;
    std::vector<int> mTrigIndex;    

    edm::Service<TFileService> fs;
    TTree *mTree,*mOutTree; 
    TFile *mInf;
    TDirectoryFile *mDir;
    //---- output TREE variables ------
    int mRun,mEvt,mNPV,mM4JIndex[2],mM2JIndex[2][2];
    float mHT,mM8J,mM4J[2],mHT4J[2],mCosThetaStar,mM4JBalance,mWeight,mPtHat;
    float mM2J[2][2],mM2JBalance[2],mDphi4J,mM2JAll[28],mM4JAll[70],mHT4JAll[70];
    float mPt[8],mEta[8],mPhi[8],mMass[8],mCHF[8],mNHF[8],mPHF[8],mELF[8],mMUF[8],mBeta[8],mBetaStar[8],mPUTag[8];
    //---- input TREE variable --------
    QCDEvent *mEvent;
};

#endif
