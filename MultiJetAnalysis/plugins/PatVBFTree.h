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
#include "TClonesArray.h"

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
    bool isBKG_;
    edm::InputTag srcJets_,srcMET_,srcRho_;
    std::string srcBtag_,srcPU_;
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    //---- output TREE variables ------
    //---- global event variables -----
    int run_,evt_,nVtx_,lumi_,nSoftTrackJets_;
    float pvx_,pvy_,pvz_,rho_,met_,metSig_,ht_,htAll_,mqq_,mbb_,dEtaqq_;
    //---- jet variables --------------
    int btagIdx_[5];
    float pt_[5],jec_[5],unc_[5],eta_[5],phi_[5],mass_[5],chf_[5],nhf_[5],phf_[5],elf_[5],muf_[5];
    float beta_[5],ptD_[5],ptMax_[5],axis_[2][5],tana_[5],ttheta_[5],btag_[5];
    TClonesArray *softTrackJetP4_;
    //---- MC variables ---------------
    int inpu_,otpu_;
    std::vector<int> *partonId_,*partonMo_,*partonSt_;
    TClonesArray *partonP4_;
};

#endif
