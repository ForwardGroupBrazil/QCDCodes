#ifndef PatVBFTree_h
#define PatVBFTree_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "TTree.h"
#include "TFile.h"

class PatVBFTree : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit PatVBFTree(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~PatVBFTree();

  private:  
    void initialize();
    //---- configurable parameters --------   
    double etaMax_,ptMin_;
    edm::InputTag srcJets_,srcMET_,srcRho_;
    std::string srcBtag_;
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    //---- output TREE variables ------
    //---- global event variables -----
    int run_,evt_,nVtx_,lumi_;
    float rho_,met_,metSig_,mqq_,mbb_,dEtaqq_;
    //---- jet variables --------------
    int btagIdx_[4];
    float pt_[4],jec_[4],eta_[4],phi_[4],mass_[4],chf_[4],nhf_[4],phf_[4],elf_[4],muf_[4];
    float beta_[4],ptD_[4],ptMax_[4],axis_[2][4],tana_[4],ttheta_[4],btag_[4];
};

#endif
