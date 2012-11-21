#ifndef PatVBFTree_h
#define PatVBFTree_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "KKousour/MultiJetAnalysis/plugins/QGLikelihoodCalculator.h"
#include "TTree.h"
#include "TH1F.h"

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
    edm::InputTag srcJets_,srcGenJets_,srcMET_,srcRho_,srcRhoQGL_;
    std::string srcBtag_,srcPU_;
    double dEtaMin_;
    std::vector<double> ptMin_;
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    //---- TRIGGER -------------------------
    triggerExpression::Data triggerCache_;
    std::vector<triggerExpression::Evaluator*> vtriggerSelector_;
    std::vector<std::string> vtriggerAlias_,vtriggerSelection_;
    TH1F *triggerPassHisto_,*triggerNamesHisto_;
    //---- output TREE variables ------
    //---- global event variables -----
    int run_,evt_,nVtx_,lumi_,nSoftTrackJets_;
    float pvx_,pvy_,pvz_,rho_,rhoQGL_,met_,metPhi_,metSig_,ht_,htAll_;
    float mqq_,mbb_,dEtaqq_,dEtabb_,dPhiqq_,dPhibb_,ptqq_,ptbb_,etaBoostqq_,etaBoostbb_,softHt_;
    float mqqEta_,mbbEta_,dEtaqqEta_;
    std::vector<bool> *triggerResult_;
    //---- jet variables --------------
    int btagIdx_[5],etaIdx_[5];
    float pt_[5],jec_[5],unc_[5],eta_[5],phi_[5],mass_[5],chf_[5],nhf_[5],phf_[5],elf_[5],muf_[5];
    float beta_[5],ptD_[5],btag_[5],qgl_[5];
    float vtxMass_[5],vtx3dL_[5],vtx3deL_[5],sumTrkPt_[5],sumTrkP_[5],sumTrkPtV_[5],leadTrkPt_[5],vtxPt_[5];
    float axis_[2][5],axis_QC_[2][5],pull_[5],pull_QC_[5],jetR_[5],jetRChg_QC_[5];
    int vtxNTrks_[5],part_[5], nChg_QC_[5],nChg_ptCut_[5],nNeutral_ptCut_[5];

    std::vector<float> *softTrackJetPt_,*softTrackJetEta_,*softTrackJetPhi_,*softTrackJetE_;
    //---- MC variables ---------------
    int npu_;
    std::vector<int> *partonId_,*partonSt_;
    std::vector<float> *partonPt_,*partonEta_,*partonPhi_,*partonE_;
    std::vector<float> *genjetPt_,*genjetEta_,*genjetPhi_,*genjetE_;
};

#endif
