#ifndef PatMultijetSearchTree_h
#define PatMultijetSearchTree_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "TTree.h"
#include "TFile.h"

class PatMultijetSearchTree : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit PatMultijetSearchTree(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~PatMultijetSearchTree();

  private:  
    std::vector<int> findIndices(const std::vector<int>& v, int rank);
    int getBin(double x, const std::vector<double>& boundaries);
    //---- configurable parameters --------   
    bool   mIsMC;
    double mEtaMax;
    double mPtMin;
    double mBetaMax;
    std::vector<double> mPtHatLumi,mPtHatBnd;
    std::string srcBeta_;
    edm::InputTag srcJets_,srcMET_;
    edm::Service<TFileService> fs;
    TTree *mOutTree; 
    //---- output TREE variables ------
    int mRun,mEvt,mNPV,mM4JIndex[2],mM2JIndex[2][2];
    float mHT,mM8J,mM4J[2],mHT4J[2],mCosThetaStar,mM4JBalance,mWeight,mPtHat,mMetSig;
    float mM2J[2][2],mM2JBalance[2],mDphi4J,mM2JAll[28],mDR2JAll[28],mM4JAll[70],mHT4JAll[70];
    float mPt[8],mEta[8],mPhi[8],mMass[8],mCHF[8],mNHF[8],mPHF[8],mELF[8],mMUF[8],mBeta[8];
};

#endif
