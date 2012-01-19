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
    std::string srcBeta_;
    edm::InputTag srcJets_,srcMET_;
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    //---- output TREE variables ------
    int run_,evt_,nVtx_,simPU_,m4JIndex_[2],m2JIndex_[2][2];
    float ht_,m8J_,m4J_[2],ht4J_[2],pt4J_[2],eta4J_[2],cosThetaStar_,m4JBalance_,metSig_;
    float m2J_[2][2],m2JBalance_[2],dPhi4J_,dR2JAll_[28];
    float pt_[8],eta_[8],phi_[8],mass_[8],chf_[8],nhf_[8],phf_[8],elf_[8],muf_[8],beta_[8];
};

#endif
