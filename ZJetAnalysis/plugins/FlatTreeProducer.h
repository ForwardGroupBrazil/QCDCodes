#ifndef FlatTreeProducer_h
#define FlatTreeProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TTree.h"
#include "TFile.h"

class  FlatTreeProducer: public edm::EDAnalyzer 
{
  public:
    typedef math::XYZTLorentzVector LorentzVector;
    explicit FlatTreeProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~FlatTreeProducer();

  private:  
    void initialize();
    //---- configurable parameters --------   
    edm::InputTag srcJets_,srcMET_,srcElectrons_,srcMuons_,srcDiElectrons_,srcDiMuons_,srcPhotons_,srcRho_;
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    double minJetPt_,maxJetEta_,minPhotonPt_;
    //---- output TREE variables ------
    int run_,evt_,lumi_,nVtx_,njets_,nelectrons_,nmuons_,nphotons_;
    float rho_;
    //---- electrons -------------------------
    std::vector<float> *elPt_,*elEta_,*elPhi_,*elE_,*elIso_;
    std::vector<int> *elId_,*elCh_;
    std::vector<LorentzVector> *elP4_;
    //---- muons -------------------------
    std::vector<float> *muPt_,*muEta_,*muPhi_,*muE_,*muIso_;
    std::vector<int> *muId_,*muCh_;
    std::vector<LorentzVector> *muP4_;
    //---- jets -------------------------
    std::vector<float> *jetPt_,*jetEta_,*jetPhi_,*jetE_,*jetBtag_,*jetRMS_,*jetBeta_,*jetChf_,*jetNhf_,*jetPhf_,*jetMuf_,*jetElf_;
    std::vector<int> *jetId_,*jetPu_;
    //---- photons -------------------------
    std::vector<float> *photonPt_,*photonEta_,*photonPhi_,*photonE_;
    //---- di-electrons -------------------------
    float eePt_,eeEta_,eePhi_,eeE_,eeM_;
    //---- di-muons -------------------------
    float mmPt_,mmEta_,mmPhi_,mmE_,mmM_;
    //---- MET -------------------------
    float met_,metPhi_;
    //---- Z/g+jet -----------------------
    float eejPt_,eejY_,eejM_,mmjPt_,mmjY_,mmjM_,gjPt_,gjY_,gjM_;
};

#endif
