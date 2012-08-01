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
    void getDoublets(const std::vector<LorentzVector>& P4, float& m2jAve, float& m2jEst, float& m2jSigma, float m2j[], int index2J[][2], int index[]);
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
    int run_,evt_,nVtx_,index_[8],index2J_[4][2],index4J_[2][4];
    float rho_,ht_,m8j_,m2j_[4],m4j_[2],ht4j_[2],pt4j_[2],eta4j_[2];
    float dPhi4j_,cosThetaStar_,m4jBalance_,metSig_,m2jAve_,m2jEst_,m2jSigma_,m4jAve_;
    float dR2jAll_[28];
    std::vector<float> *pt_,*eta_,*phi_,*mass_,*chf_,*nhf_,*phf_,*elf_,*muf_,*beta_,*jec_,*unc_;
    //---- MC variables ---------------
    int simPU_;
    std::vector<int> *partonId_,*partonSt_;
    std::vector<float> *partonPt_,*partonEta_,*partonPhi_,*partonE_;
    float m2jAveGEN_,m2jEstGEN_,m4jAveGEN_,m2jAveParton_,m2jEstParton_,m4jAveParton_;
};

#endif
