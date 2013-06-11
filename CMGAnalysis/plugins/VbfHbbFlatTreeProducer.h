#ifndef VbfHbbFlatTreeProducer_h
#define VbfHbbFlatTreeProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "TTree.h"
#include "TH1F.h"
#include "KKousour/CMGAnalysis/plugins/QGLCalculator.h"
#include "KKousour/CMGAnalysis/plugins/JetRegressor.h"

class VbfHbbFlatTreeProducer : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit VbfHbbFlatTreeProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~VbfHbbFlatTreeProducer();

  private:  
    void initialize();
    void order(std::vector<float> const& v, std::vector<int>* idx);
    //---- configurable parameters --------   
    edm::InputTag srcJets_,srcGenJets_,srcMET_,srcRho_,srcGenParticles_;
    std::string srcBtag_,srcPU_;
    double dEtaMin_,shiftJES_;
    double ptMin_;
    std::vector<double> ptPreselMin_;
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    //---- TRIGGER -------------------------
    triggerExpression::Data triggerCache_;
    std::vector<triggerExpression::Evaluator*> vtriggerSelector_;
    std::vector<std::string> vtriggerAlias_,vtriggerSelection_;
    TH1F *triggerPassHisto_,*triggerNamesHisto_;
    //---- output TREE variables ------
    //---- global event variables -----
    int run_,evt_,nVtx_,lumi_,nSoftTrackJets_,nJets_,nBJets_;
    int b1_,b2_,q1_,q2_;
    float pvx_,pvy_,pvz_,rho_,met_,metPhi_,metSig_,ht_,htAll_,sphericity_,aplanarity_,cosTheta_;
    float mqq_,mbb_,mbbReg_,mbbRegPart_,dEtaqq_,dEtabb_,dPhiqq_,dPhibb_,ptqq_,ptbb_,etaBoostqq_,etaBoostbb_,softHt_,dEtaMax_;
    std::vector<bool> *triggerResult_;
    //---- jet variables --------------
    std::vector<bool>  *puIdL_,*puIdM_,*puIdT_,*idL_,*idM_,*idT_;
    std::vector<int>   *btagIdx_,*etaIdx_;
    std::vector<int>   *vtxNTrks_,*part_,*nChg_QC_,*nChg_ptCut_,*nNeutral_ptCut_;
    std::vector<float> *pt_,*jec_,*reg_,*regPart_,*unc_,*eta_,*phi_,*mass_,*chf_,*nhf_,*phf_,*elf_,*muf_,*jetMetPhi_;
    std::vector<float> *beta_,*ptD_,*ptD_QC_,*btag_,*puMva_,*qgl_;
    std::vector<float> *vtxPt_,*vtx3dL_,*vtx3deL_;
    std::vector<float> *axisMinor_,*axisMajor_,*axisMinor_QC_,*axisMajor_QC_,*pull_,*pull_QC_,*jetR_,*jetRChg_QC_;

    std::vector<float> *softTrackJetPt_,*softTrackJetEta_,*softTrackJetPhi_,*softTrackJetE_;
    //---- MC variables ---------------
    int npu_;
    std::vector<int>   *partonId_,*partonSt_,*partonMatchIdx_;
    std::vector<float> *partonPt_,*partonEta_,*partonPhi_,*partonE_,*partonMatchDR_;
    std::vector<float> *genjetPt_,*genjetEta_,*genjetPhi_,*genjetE_;
    //---- QGL tagger -----------------
    QGLCalculator *qglCalc_;
    //---- Jet regressor -----------------
    JetRegressor *jetReg_,*jetRegPart_;
};

#endif
