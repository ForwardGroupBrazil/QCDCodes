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
    void initialize();
    std::vector<int> findIndices(const std::vector<int>& v, int rank);
    //---- configurable parameters --------   
    double etaMax_,ptMin_,betaMax_;
    std::string srcBeta_,srcPU_;
    edm::InputTag srcJets_,srcMET_,srcRho_;
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    //---- output TREE variables ------
    int run_,evt_,nVtx_,simPU_,index2J_[4][2],index4J_[2][4];
    float rho_,ht_,m8j_,m2j_[4],m4j_[2],ht4j_[2],pt4j_[2],eta4j_[2],dPhi4j_,cosThetaStar_,m4jBalance_,metSig_,m2jAve_,m2jSigma_,m4jAve_;
    float dR2jAll_[28];
    float pt_[8],eta_[8],phi_[8],mass_[8],chf_[8],nhf_[8],phf_[8],elf_[8],muf_[8],beta_[8];
};

#endif
