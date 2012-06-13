#ifndef SimpleTreeProducer_h
#define SimpleTreeProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "TTree.h"
#include "TFile.h"

class  SimpleTreeProducer: public edm::EDAnalyzer 
{
  public:
    explicit SimpleTreeProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~SimpleTreeProducer();

  private:  
    void initialize();
    
    //---- configurable parameters --------   
    bool isMC_;
    std::string srcPU_;
    edm::InputTag srcJets_,srcMET_,srcRho_,srcGenJets_;
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    //---- output TREE variables ------
    int run_,evt_,nVtx_;
    //---- lepton -------------------------
    float recoLeptonPt_,recoLeptonEta_,recoLeptonPhi_,recoLeptonE_;
    float genLeptonPt_,genLeptonEta_,genLeptonPhi_,genLeptonE_;
    //---- BJet -------------------------
    float recoBJetPt_,recoBJetEta_,recoBJetPhi_,recoBJetE_;
    float genBJetPt_,genBJetEta_,genBJetPhi_,genBJetE_;
    //---- QJet -------------------------
    float recoQJetPt_,recoQJetEta_,recoQJetPhi_,recoQJetE_;
    float genQJetPt_,genQJetEta_,genQJetPhi_,genQJetE_;
    //---- MET -------------------------
    float recoMet_,recoMetPhi_;
    float genMet_,genMetPhi_;
    //---- W ----------------------------
    float partonWPt_,partonWEta_,partonWPhi_,partonWE_;
    //---- TOP --------------------------
    float partonTopPt_,partonTopEta_,partonTopPhi_,partonTopE_;
    
};

#endif
