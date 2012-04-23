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
    void getDoublets(const std::vector<LorentzVector>& P4, float& m2jAve, float& m2jSigma, float m2j[], int index2J[][2]);
    void getQuartets(const std::vector<LorentzVector>& P4, LorentzVector quartetP4[], float& m4jAve, float& m4jBalance, float m4j[], int index4J[][4], int index2J[][2]);
    std::vector<int> findIndices(const std::vector<int>& v, int rank);
    //---- configurable parameters --------   
    bool isMC_;
    double etaMax_,ptMin_,betaMin_;
    std::string srcPU_;
    edm::InputTag srcJets_,srcMET_,srcRho_,srcGenJets_;
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    //---- output TREE variables ------
    int run_,evt_,nVtx_,index2J_[4][2],index4J_[2][4];
    float rho_,ht_,m8j_,m2j_[4],m4j_[2],ht4j_[2],pt4j_[2],eta4j_[2],dPhi4j_,cosThetaStar_,m4jBalance_,metSig_,m2jAve_,m2jSigma_,m4jAve_;
    float dR2jAll_[28];
    float pt_[8],eta_[8],phi_[8],mass_[8],chf_[8],nhf_[8],phf_[8],elf_[8],muf_[8],beta_[8],jec_[8],unc_[8];
    //---- MC variables ---------------
    int simPU_;
    std::vector<int> *partonId_,*partonSt_;
    std::vector<float> *partonPt_,*partonEta_,*partonPhi_,*partonE_;
    float m2jAveGEN_,m4jAveGEN_,m2jAveParton_,m4jAveParton_;
};

#endif
