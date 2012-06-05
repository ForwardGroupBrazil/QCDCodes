#ifndef LowerMultiplicityTree_h
#define LowerMultiplicityTree_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "TTree.h"
#include "TFile.h"

class LowerMultiplicityTree : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit LowerMultiplicityTree(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~LowerMultiplicityTree();

  private:  
    void initialize();
    
    //---- configurable parameters --------   
    bool isMC_;
    double etaMax_,ptMin_,betaMin_;
    std::string srcPU_;
    edm::InputTag srcJets_,srcMET_,srcRho_,srcGenJets_;
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    //---- output TREE variables ------
    int run_,evt_,nVtx_;
    float ht_;
    float pt_[8],eta_[8],phi_[8];
    int njet_;
    
};

#endif
