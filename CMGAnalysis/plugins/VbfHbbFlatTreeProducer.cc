#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "TMath.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "KKousour/CMGAnalysis/plugins/VbfHbbFlatTreeProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "AnalysisDataFormats/CMGTools/interface/BaseMET.h"
#include "CMGTools/External/interface/PileupJetIdentifier.h"

using namespace std;
using namespace reco;

VbfHbbFlatTreeProducer::VbfHbbFlatTreeProducer(edm::ParameterSet const& cfg) 
{
  srcJets_            = cfg.getParameter<edm::InputTag>             ("jets");
  srcMET_             = cfg.getParameter<edm::InputTag>             ("met");
  srcRho_             = cfg.getParameter<edm::InputTag>             ("rho"); 
  srcBtag_            = cfg.getParameter<std::string>               ("btagger");
  shiftJES_           = cfg.getParameter<double>                    ("shiftJES");
  ptMin_              = cfg.getParameter<double>                    ("ptMin");
  dEtaMin_            = cfg.getParameter<double>                    ("dEtaMin");
  srcPU_              = cfg.getUntrackedParameter<std::string>      ("pu","");
  srcGenJets_         = cfg.getUntrackedParameter<edm::InputTag>    ("genjets",edm::InputTag(""));
  srcGenParticles_    = cfg.getUntrackedParameter<edm::InputTag>    ("genparticles",edm::InputTag(""));

  triggerCache_       = triggerExpression::Data(cfg.getParameterSet("triggerConfiguration"));
  vtriggerAlias_      = cfg.getParameter<std::vector<std::string> > ("triggerAlias");
  vtriggerSelection_  = cfg.getParameter<std::vector<std::string> > ("triggerSelection");

  if (vtriggerAlias_.size() != vtriggerSelection_.size()) {
    cout<<"ERROR: the number of trigger aliases does not match the number of trigger names !!!"<<endl;
    return;
  }

  for(unsigned i=0;i<vtriggerSelection_.size();i++) {
    vtriggerSelector_.push_back(triggerExpression::parse(vtriggerSelection_[i]));
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void VbfHbbFlatTreeProducer::beginJob() 
{
  //--- book the trigger histograms ---------
  triggerNamesHisto_ = fs_->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  triggerNamesHisto_->SetBit(TH1::kCanRebin);
  for(unsigned i=0;i<vtriggerSelection_.size();i++) {
    triggerNamesHisto_->Fill(vtriggerSelection_[i].c_str(),1);
  }
  triggerPassHisto_ = fs_->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  triggerPassHisto_->SetBit(TH1::kCanRebin);
  //--- book the tree -----------------------
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"                ,&run_               ,"run_/I");
  outTree_->Branch("evtNo"                ,&evt_               ,"evt_/I");
  outTree_->Branch("lumi"                 ,&lumi_              ,"lumi_/I");
  outTree_->Branch("nvtx"                 ,&nVtx_              ,"nVtx_/I");
  outTree_->Branch("nSoftTrackJets"       ,&nSoftTrackJets_    ,"nSoftTrackJets_/I");
  outTree_->Branch("nJets"                ,&nJets_             ,"nJets_/I");
  outTree_->Branch("nBJets"               ,&nBJets_            ,"nBJets_/I");
  outTree_->Branch("b1"                   ,&b1_                ,"b1_/I");
  outTree_->Branch("b2"                   ,&b2_                ,"b2_/I");
  outTree_->Branch("q1"                   ,&q1_                ,"q1_/I");
  outTree_->Branch("q2"                   ,&q2_                ,"q2_/I"); 
  outTree_->Branch("softHt"               ,&softHt_            ,"softHt_/F");
  outTree_->Branch("pvx"                  ,&pvx_               ,"pvx_/F");
  outTree_->Branch("pvy"                  ,&pvy_               ,"pvy_/F");
  outTree_->Branch("pvz"                  ,&pvz_               ,"pvz_/F");
  outTree_->Branch("rho"                  ,&rho_               ,"rho_/F");
  outTree_->Branch("ht"                   ,&ht_                ,"ht_/F");
  outTree_->Branch("pvz"                  ,&pvz_               ,"pvz_/F");
  outTree_->Branch("met"                  ,&met_               ,"met_/F");
  outTree_->Branch("metPhi"               ,&metPhi_            ,"metPhi_/F");
  outTree_->Branch("metSig"               ,&metSig_            ,"metSig_/F");
  outTree_->Branch("sphericity"           ,&sphericity_        ,"sphericity_/F");
  outTree_->Branch("aplanarity"           ,&aplanarity_        ,"aplanarity_/F");
  outTree_->Branch("mbb"                  ,&mbb_               ,"mbb_/F");
  outTree_->Branch("mbbNew"               ,&mbbNew_            ,"mbbNew_/F");
  outTree_->Branch("dEtaMax"              ,&dEtaMax_           ,"dEtaMax_/F");
  outTree_->Branch("dEtaqq"               ,&dEtaqq_            ,"dEtaqq_/F");
  //------------------------------------------------------------------
  btagIdx_        = new std::vector<int>;
  etaIdx_         = new std::vector<int>;
  puIdL_          = new std::vector<bool>;
  puIdM_          = new std::vector<bool>;
  puIdT_          = new std::vector<bool>;
  idL_            = new std::vector<bool>;
  idM_            = new std::vector<bool>;
  idT_            = new std::vector<bool>;
  pt_             = new std::vector<float>;
  btag_           = new std::vector<float>;
  puMva_          = new std::vector<float>;
  jec_            = new std::vector<float>;
  unc_            = new std::vector<float>;
  beta_           = new std::vector<float>;
  eta_            = new std::vector<float>;
  phi_            = new std::vector<float>;
  mass_           = new std::vector<float>;
  chf_            = new std::vector<float>;
  nhf_            = new std::vector<float>;
  phf_            = new std::vector<float>;
  muf_            = new std::vector<float>;
  elf_            = new std::vector<float>;
  ptD_            = new std::vector<float>;
  ptD_QC_         = new std::vector<float>;
  axisMinor_      = new std::vector<float>;
  axisMajor_      = new std::vector<float>;
  axisMinor_QC_   = new std::vector<float>;
  axisMajor_QC_   = new std::vector<float>;
  pull_           = new std::vector<float>;
  pull_QC_        = new std::vector<float>;
  jetR_           = new std::vector<float>;
  jetRChg_QC_     = new std::vector<float>;
  vtx3dL_         = new std::vector<float>;
  vtx3deL_        = new std::vector<float>;
  vtxPt_          = new std::vector<float>;
  part_           = new std::vector<int>;
  nChg_QC_        = new std::vector<int>;
  nChg_ptCut_     = new std::vector<int>;
  nNeutral_ptCut_ = new std::vector<int>;
  outTree_->Branch("btagIdx"              ,"vector<int>"       ,&btagIdx_);
  outTree_->Branch("etaIdx"               ,"vector<int>"       ,&etaIdx_);
  outTree_->Branch("jetPuIdL"             ,"vector<bool>"      ,&puIdL_);
  outTree_->Branch("jetPuIdM"             ,"vector<bool>"      ,&puIdM_);
  outTree_->Branch("jetPuIdT"             ,"vector<bool>"      ,&puIdT_);
  outTree_->Branch("jetIdL"               ,"vector<bool>"      ,&idL_);
  outTree_->Branch("jetIdM"               ,"vector<bool>"      ,&idM_);
  outTree_->Branch("jetIdT"               ,"vector<bool>"      ,&idT_); 
  outTree_->Branch("jetPt"                ,"vector<float>"     ,&pt_);
  outTree_->Branch("jetBtag"              ,"vector<float>"     ,&btag_);
  outTree_->Branch("jetPuMva"             ,"vector<float>"     ,&puMva_);
  outTree_->Branch("jetJec"               ,"vector<float>"     ,&jec_);
  outTree_->Branch("jetUnc"               ,"vector<float>"     ,&unc_);
  outTree_->Branch("jetBeta"              ,"vector<float>"     ,&beta_);
  outTree_->Branch("jetEta"               ,"vector<float>"     ,&eta_);
  outTree_->Branch("jetPhi"               ,"vector<float>"     ,&phi_);
  outTree_->Branch("jetMass"              ,"vector<float>"     ,&mass_);
  outTree_->Branch("jetChf"               ,"vector<float>"     ,&chf_);
  outTree_->Branch("jetNhf"               ,"vector<float>"     ,&nhf_);
  outTree_->Branch("jetPhf"               ,"vector<float>"     ,&phf_);
  outTree_->Branch("jetMuf"               ,"vector<float>"     ,&muf_);
  outTree_->Branch("jetElf"               ,"vector<float>"     ,&elf_);
  outTree_->Branch("jetPtD"               ,"vector<float>"     ,&ptD_);
  outTree_->Branch("jetPtD_QC"            ,"vector<float>"     ,&ptD_QC_);
  //------------------------------------------------------------------- 
  outTree_->Branch("jetAxisMinor"         ,"vector<float>"     ,&axisMinor_);
  outTree_->Branch("jetAxisMajor"         ,"vector<float>"     ,&axisMajor_);
  outTree_->Branch("jetAxisMinor_QC"      ,"vector<float>"     ,&axisMinor_QC_);
  outTree_->Branch("jetAxisMajor_QC"      ,"vector<float>"     ,&axisMajor_QC_);
  outTree_->Branch("jetPull"              ,"vector<float>"     ,&pull_);
  outTree_->Branch("jetPull_QC"           ,"vector<float>"     ,&pull_QC_);
  outTree_->Branch("jetR"                 ,"vector<float>"     ,&jetR_);
  outTree_->Branch("jetRChg_QC"           ,"vector<float>"     ,&jetRChg_QC_);
  outTree_->Branch("jetVtx3dL"            ,"vector<float>"     ,&vtx3dL_);
  outTree_->Branch("jetVtx3deL"           ,"vector<float>"     ,&vtx3deL_);
  outTree_->Branch("jetVtxPt"             ,"vector<float>"     ,&vtxPt_);
  outTree_->Branch("jetPart"              ,"vector<int>"       ,&part_); 
  outTree_->Branch("jetnChg_QC"           ,"vector<int>"       ,&nChg_QC_);
  outTree_->Branch("jetnChg_ptCut"        ,"vector<int>"       ,&nChg_ptCut_);
  outTree_->Branch("jetnNeutral_ptCut"    ,"vector<int>"       ,&nNeutral_ptCut_);
  //------------------------------------------------------------------
  triggerResult_ = new std::vector<bool>;
  outTree_->Branch("triggerResult","vector<bool>",&triggerResult_);
  //------------------------------------------------------------------
  softTrackJetPt_  = new std::vector<float>;
  softTrackJetEta_ = new std::vector<float>;
  softTrackJetPhi_ = new std::vector<float>;
  softTrackJetE_   = new std::vector<float>;
  outTree_->Branch("softTrackJetPt" ,"vector<float>" ,&softTrackJetPt_);
  outTree_->Branch("softTrackJetEta","vector<float>" ,&softTrackJetEta_);
  outTree_->Branch("softTrackJetPhi","vector<float>" ,&softTrackJetPhi_);
  outTree_->Branch("softTrackJetE"  ,"vector<float>" ,&softTrackJetE_);
  //------------------- MC ---------------------------------
  outTree_->Branch("npu"           ,&npu_          ,"npu_/I");
  partonId_       = new std::vector<int>;
  partonSt_       = new std::vector<int>;
  partonMatchIdx_ = new std::vector<int>;
  partonMatchDR_  = new std::vector<float>;
  partonPt_       = new std::vector<float>;
  partonEta_      = new std::vector<float>;
  partonPhi_      = new std::vector<float>;
  partonE_        = new std::vector<float>;
  genjetPt_       = new std::vector<float>;
  genjetEta_      = new std::vector<float>;
  genjetPhi_      = new std::vector<float>;
  genjetE_        = new std::vector<float>;
  outTree_->Branch("partonId"       ,"vector<int>"   ,&partonId_);
  outTree_->Branch("partonSt"       ,"vector<int>"   ,&partonSt_);
  outTree_->Branch("partonMatchIdx" ,"vector<int>"   ,&partonMatchIdx_);
  outTree_->Branch("partonMatchDR"  ,"vector<float>" ,&partonMatchDR_);
  outTree_->Branch("partonPt"       ,"vector<float>" ,&partonPt_);
  outTree_->Branch("partonEta"      ,"vector<float>" ,&partonEta_);
  outTree_->Branch("partonPhi"      ,"vector<float>" ,&partonPhi_);
  outTree_->Branch("partonE"        ,"vector<float>" ,&partonE_);
  outTree_->Branch("genjetPt"       ,"vector<float>" ,&genjetPt_);
  outTree_->Branch("genjetEta"      ,"vector<float>" ,&genjetEta_);
  outTree_->Branch("genjetPhi"      ,"vector<float>" ,&genjetPhi_);
  outTree_->Branch("genjetE"        ,"vector<float>" ,&genjetE_);

  outTree_->Branch("partonDRbb"     ,&partonDRbb_    ,"partonDRbb_/F");

  cout<<"Begin job finished"<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void VbfHbbFlatTreeProducer::endJob() 
{
  delete partonSt_;
  delete partonId_;
  delete partonMatchIdx_;
  delete partonMatchDR_;
  delete partonPt_;
  delete partonEta_;
  delete partonPhi_;
  delete partonE_;
  delete genjetPt_;
  delete genjetEta_;
  delete genjetPhi_;
  delete genjetE_;
  delete softTrackJetPt_;
  delete softTrackJetEta_;
  delete softTrackJetPhi_;
  delete softTrackJetE_;
  delete triggerResult_;
  delete btagIdx_;
  delete etaIdx_;
  delete puIdL_;
  delete puIdM_;
  delete puIdT_;
  delete idL_;
  delete idM_;
  delete idT_;
  delete pt_;
  delete btag_;
  delete puMva_;
  delete jec_;
  delete unc_;
  delete beta_;
  delete eta_;
  delete phi_;
  delete mass_;
  delete chf_;
  delete nhf_;
  delete phf_;
  delete muf_;
  delete elf_;
  delete ptD_;
  delete ptD_QC_;
  delete axisMinor_;
  delete axisMajor_;
  delete axisMinor_QC_;
  delete axisMajor_QC_;
  delete pull_;
  delete pull_QC_;
  delete jetR_;
  delete jetRChg_QC_;
  delete vtx3dL_;
  delete vtx3deL_;
  delete vtxPt_;
  delete part_;
  delete nChg_QC_;
  delete nChg_ptCut_;
  delete nNeutral_ptCut_;
  
  for(unsigned i=0;i<vtriggerSelector_.size();i++) {
    delete vtriggerSelector_[i];
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void VbfHbbFlatTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  edm::Handle<edm::View<cmg::PFJet> > jets;
  iEvent.getByLabel(srcJets_,jets);
  edm::View<cmg::PFJet> cmg_jets = *jets;

  edm::Handle<reco::TrackJetCollection> softTrackJets;
  iEvent.getByLabel("ak5SoftTrackJetsForVbfHbb",softTrackJets);

  edm::Handle<edm::View<cmg::BaseMET> >  met;
  iEvent.getByLabel(srcMET_,met);

  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("offlinePrimaryVertices",recVtxs);

  edm::Handle<GenEventInfoProduct> hEventInfo;
  edm::Handle<GenParticleCollection> genParticles;
  
  edm::Handle<std::vector<cmg::PhysicsObjectWithPtr<edm::Ptr<reco::GenJet> > > > genjets;
  //edm::Handle<reco::GenJetCollection> genjets;
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;

  initialize();
 
  //-------------- Trigger Info -----------------------------------
  triggerPassHisto_->Fill("totalEvents",1);
  if (triggerCache_.setEvent(iEvent,iSetup)) {
    for(unsigned itrig=0;itrig<vtriggerSelector_.size();itrig++) {
      bool result(false);
      if (vtriggerSelector_[itrig]) {
        if (triggerCache_.configurationUpdated()) {
          vtriggerSelector_[itrig]->init(triggerCache_);
        }
        result = (*(vtriggerSelector_[itrig]))(triggerCache_);
      }
      if (result) {
        triggerPassHisto_->Fill(vtriggerAlias_[itrig].c_str(),1);
      }
      triggerResult_->push_back(result);
    }
  }
     
  //----- at least one good vertex -----------
  bool cut_vtx = (recVtxs->size() > 0);
  
  if (cut_vtx) {
    //----- soft track jets ----------------------
    float softHt(0.0);
    int nsoft(0); 
    for(int it=0;it<int(softTrackJets->size());it++) {
      //--- make sure that the track jets do not match with the 4 leading PFJets ----
      int njets(0);
      bool matched(false);
      for(edm::View<cmg::PFJet>::const_iterator ijet = cmg_jets.begin();(ijet != cmg_jets.end() && njets<4); ++ijet) {
        double dR = deltaR(ijet->p4(),(*softTrackJets)[it].p4());
        if (dR < 0.5) {
          matched = true;
          break;
        }
        njets++;
      }
      if (matched) continue;
      softHt += (*softTrackJets)[it].pt();
      if (nsoft < 3) {
        softTrackJetPt_ ->push_back((*softTrackJets)[it].pt());
        softTrackJetEta_->push_back((*softTrackJets)[it].eta());
        softTrackJetPhi_->push_back((*softTrackJets)[it].phi());
        softTrackJetE_  ->push_back((*softTrackJets)[it].energy());
      } 
      nsoft++;
    }
    nSoftTrackJets_ = nsoft;
    //----- PF jets ------------------------------
    nJets_ = 0;
    nBJets_ = 0;
    float ht(0.0);
    float sumP2(0.0),sumPxx(0.0),sumPxy(0.0),sumPxz(0.0),sumPyy(0.0),sumPyz(0.0),sumPzz(0.0);
    vector<float> vaJES;
    vector<TLorentzVector> vP4;
    for(edm::View<cmg::PFJet>::const_iterator ijet = cmg_jets.begin();ijet != cmg_jets.end(); ++ijet) {  
      float chf = ijet->component(1).fraction();
      float nhf = ijet->component(5).fraction();
      float phf = ijet->component(4).fraction();
      float muf = ijet->component(3).fraction();
      float elf = ijet->component(2).fraction();
      int chm = ijet->component(1).number();
      int npr = ijet->nConstituents();
      float eta = fabs(ijet->eta());
      float pt = ijet->pt();
      bool puIdL = ijet->passPuJetId("full",PileupJetIdentifier::kLoose);
      bool puIdM = ijet->passPuJetId("full",PileupJetIdentifier::kMedium);
      bool puIdT = ijet->passPuJetId("full",PileupJetIdentifier::kTight);
      bool idL = (npr>1 && phf<0.99 && nhf<0.99);
      bool idM = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || eta>2.4));
      bool idT = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.7 && muf<0.7 && chf>0 && chm>0) || eta>2.4));
      if (idL && puIdL && (pt > ptMin_)) {
        vP4.push_back(TLorentzVector(ijet->px(),ijet->py(),ijet->pz(),ijet->energy()));
        sumP2  += ijet->p() * ijet->p();
        sumPxx += ijet->px() * ijet->px();
        sumPxy += ijet->px() * ijet->py();
        sumPxz += ijet->px() * ijet->pz();
        sumPyy += ijet->py() * ijet->py();
        sumPyz += ijet->py() * ijet->pz();
        sumPzz += ijet->pz() * ijet->pz();   
        double aJES = 1.+shiftJES_*ijet->uncOnFourVectorScale();
        vaJES.push_back(aJES);
        chf_           ->push_back(chf);
        nhf_           ->push_back(nhf);
        phf_           ->push_back(phf);
        elf_           ->push_back(elf);
        muf_           ->push_back(muf);
        unc_           ->push_back(ijet->uncOnFourVectorScale());
        jec_           ->push_back(aJES/ijet->rawFactor());
        pt_            ->push_back(aJES*pt);
        phi_           ->push_back(ijet->phi());
        eta_           ->push_back(ijet->eta());
        mass_          ->push_back(aJES*ijet->mass());
        btag_          ->push_back(ijet->bDiscriminator(srcBtag_.c_str()));
        beta_          ->push_back(ijet->beta());
        puMva_         ->push_back(ijet->puMva("full"));
        puIdL_         ->push_back(puIdL);
        puIdM_         ->push_back(puIdM);
        puIdT_         ->push_back(puIdT);
        idL_           ->push_back(idL);
        idM_           ->push_back(idM);
        idT_           ->push_back(idT);
        ptD_           ->push_back(ijet->ptd());
        ptD_QC_        ->push_back(ijet->ptdQC());
        vtx3dL_        ->push_back(ijet->vtx3dL());
        vtx3deL_       ->push_back(ijet->vtx3deL());
        vtxPt_         ->push_back(ijet->vtxPt());
        part_          ->push_back(npr);
        axisMajor_     ->push_back(ijet->axisMajor());
        axisMinor_     ->push_back(ijet->axisMinor());
        axisMajor_QC_  ->push_back(ijet->axisMajorQC());
        axisMinor_QC_  ->push_back(ijet->axisMinorQC());
        pull_          ->push_back(ijet->pull());
        pull_QC_       ->push_back(ijet->pullQC());
        jetR_          ->push_back(ijet->fmax());
        jetRChg_QC_    ->push_back(ijet->fmaxCharged());
        nChg_ptCut_    ->push_back(ijet->nChargedPtCut());
        nNeutral_ptCut_->push_back(ijet->nNeutralPtCut());
        nChg_QC_       ->push_back(ijet->nChargedQC());
        ht += aJES*pt;
        nJets_++;
        if (ijet->bDiscriminator(srcBtag_.c_str()) > 0.679) {
          nBJets_++; 
        }
      }
    }// jet loop
    //---- compute sphericity -----------------------
    float Txx = sumPxx/sumP2;
    float Tyy = sumPyy/sumP2;
    float Tzz = sumPzz/sumP2;  
    float Txy = sumPxy/sumP2; 
    float Txz = sumPxz/sumP2; 
    float Tyz = sumPyz/sumP2;
    TMatrixDSym T(3);
    T(0,0) = Txx;
    T(0,1) = Txy;
    T(0,2) = Txz;
    T(1,0) = Txy;
    T(1,1) = Tyy;
    T(1,2) = Tyz;
    T(2,0) = Txz;
    T(2,1) = Tyz;
    T(2,2) = Tzz;
    TMatrixDSymEigen TEigen(T);
    TVectorD eigenValues(TEigen.GetEigenValues());
    sphericity_ = 1.5*(eigenValues(1)+eigenValues(2));
    aplanarity_ = 1.5*eigenValues(2); 
    //----- order jets according to the btag value -----
    std::vector<float> tmp_btag;
    for(int i=0;i<nJets_;i++) {
      btagIdx_->push_back(i);
      tmp_btag.push_back((*btag_)[i]); 
    }
    for(int i=0;i<nJets_;i++) {
      for(int j=i+1;j<nJets_;j++) {
        if (tmp_btag[j] > tmp_btag[i]) {
          int k = (*btagIdx_)[i];
          (*btagIdx_)[i] = (*btagIdx_)[j];
          (*btagIdx_)[j] = k;
          float tmp = tmp_btag[i]; 
          tmp_btag[i] = tmp_btag[j];
          tmp_btag[j] = tmp;
        }
      } 
    }
    //----- order jets according to the eta value -----
    std::vector<float> tmp_eta;
    for(int i=0;i<nJets_;i++) {
      etaIdx_->push_back(i);
      tmp_eta.push_back((*eta_)[i]); 
    }
    for(int i=0;i<nJets_;i++) {
      for(int j=i+1;j<nJets_;j++) {
        if (tmp_eta[j] > tmp_eta[i]) {
          int k = (*etaIdx_)[i];
          (*etaIdx_)[i] = (*etaIdx_)[j];
          (*etaIdx_)[j] = k;
          float tmp = tmp_eta[i]; 
          tmp_eta[i] = tmp_eta[j];
          tmp_eta[j] = tmp;
        }
      } 
    }
    // ------- MC -----------------------------------------
    if (!iEvent.isRealData()) {
      //---------- pu -----------------------
      if (srcPU_ != "") {
        iEvent.getByLabel(srcPU_,PupInfo);
        std::vector<PileupSummaryInfo>::const_iterator PUI;
        for(PUI = PupInfo->begin(); PUI != PupInfo->end(); ++PUI) {
          if (PUI->getBunchCrossing() == 0) {
            npu_ = PUI->getTrueNumInteractions();
          }
        }
      }
      //---------- genjets ------------------
      iEvent.getByLabel(srcGenJets_, genjets);
      //for(reco::GenJetCollection::const_iterator igen = genjets->begin();igen != genjets->end(); ++igen) {
      for(std::vector<cmg::PhysicsObjectWithPtr<edm::Ptr<reco::GenJet> > >::const_iterator igen = genjets->begin();igen != genjets->end(); ++igen) {
        genjetPt_ ->push_back(igen->pt());
        genjetEta_->push_back(igen->eta());
        genjetPhi_->push_back(igen->phi());
        genjetE_  ->push_back(igen->energy()); 
      }// genjet loop
      //---------- partons ------------------
      iEvent.getByLabel(srcGenParticles_, genParticles);
      for(unsigned ip = 0; ip < genParticles->size(); ++ ip) {
        const GenParticle &p = (*genParticles)[ip];
        /*
        if (p.pdgId() == 25) {
          partonId_ ->push_back(p.pdgId());
	  partonSt_ ->push_back(p.status()); 
          partonPt_ ->push_back(p.pt());
          partonEta_->push_back(p.eta());
          partonPhi_->push_back(p.phi());
          partonE_  ->push_back(p.energy());
          partonDRbb_ = deltaR(p.daughter(0)->p4(),p.daughter(1)->p4());
        }
        */
        if (p.status() != 3) continue;
        int ndst3(0);
	for(unsigned k = 0; k < p.numberOfDaughters(); k++) {
	  if (p.daughter(k)->status() == 3) {
            ndst3++; 
          }
	}
        if (ndst3 > 0) continue;
        partonId_ ->push_back(p.pdgId());
	partonSt_ ->push_back(p.status()); 
        partonPt_ ->push_back(p.pt());
        partonEta_->push_back(p.eta());
        partonPhi_->push_back(p.phi());
        partonE_  ->push_back(p.energy()); 
        //----- match partons with jets ------------
        float dRmin(1000);
        float eta1 = p.eta();
        float phi1 = p.phi();  
        int imatch(0);
        for(int j=0;j<nJets_;j++) {
          float eta2 = (*eta_)[j];
          float phi2 = (*phi_)[j];
          float de = std::abs(eta1-eta2);
          float dp = std::abs(phi1-phi2); 
          if (dp > 3.14159) dp -= (2*3.14159);  
          float dR = sqrt(de*de + dp*dp);
          if (dR < dRmin) {
            imatch = j;
            dRmin = dR;
          }
        }   
        partonMatchIdx_->push_back(imatch);
        partonMatchDR_->push_back(dRmin);
      }// parton loop
    }// if MC   
    ht_     = ht;
    softHt_ = softHt;
    rho_    = *rho;
    pvx_    = (*recVtxs)[0].x();
    pvy_    = (*recVtxs)[0].y();
    pvz_    = (*recVtxs)[0].z();
    nVtx_   = recVtxs->size();
    met_    = (*met)[0].et();
    metPhi_ = (*met)[0].phi();
    metSig_ = (*met)[0].et()/(*met)[0].sumEt();
    run_    = iEvent.id().run();
    evt_    = iEvent.id().event();
    lumi_   = iEvent.id().luminosityBlock();
    if (nJets_ > 3) {
      b1_ = (*btagIdx_)[0];
      b2_ = (*btagIdx_)[1];
      q1_ = (*btagIdx_)[2];
      q2_ = (*btagIdx_)[3]; 
      float btag21 = (*btag_)[(*btagIdx_)[2]]/(*btag_)[(*btagIdx_)[1]];
      if (btag21 > 0.5) {
        if ((*etaIdx_)[(*btagIdx_)[1]] == 0 || (*etaIdx_)[(*btagIdx_)[1]] > 2) {
          b2_ = (*btagIdx_)[2];
          q1_ = (*btagIdx_)[1];
          q2_ = (*btagIdx_)[3];
        }
      } 
      int ib0 = (*btagIdx_)[0];
      int ib1 = (*btagIdx_)[1];
      mbb_    = (vaJES[ib0]*vP4[ib0]+vaJES[ib1]*vP4[ib1]).M();
      mbbNew_ = (vaJES[b1_]*vP4[b1_]+vaJES[b2_]*vP4[b2_]).M(); 
      cout<<btag21<<" "<<mbb_<<" "<<mbbNew_<<" "<<ib0<<" "<<ib1<<" "<<b1_<<" "<<b2_<<endl;
      int ie0 = (*etaIdx_)[0];
      int ie1 = (*etaIdx_)[nJets_-1];
      dEtaMax_ = fabs((*eta_)[ie0]-(*eta_)[ie1]);
      dEtaqq_  = fabs((*eta_)[q1_]-(*eta_)[q2_]);
      outTree_->Fill();
    }
  }// if vtx
}
//////////////////////////////////////////////////////////////////////////////////////////
void VbfHbbFlatTreeProducer::initialize()
{
  run_            = -999;
  evt_            = -999;
  lumi_           = -999;
  nVtx_           = -999;
  nSoftTrackJets_ = -999;
  nJets_          = -999;
  nBJets_         = -999;
  rho_            = -999;
  met_            = -999;
  metPhi_         = -999;
  metSig_         = -999;
  ht_             = -999;
  softHt_         = -999;
  pvx_            = -999;
  pvy_            = -999;
  pvz_            = -999;
  sphericity_     = -999;
  aplanarity_     = -999;
  mbb_            = -999;
  mbbNew_         = -999;
  dEtaMax_        = -999;
  dEtaqq_         = -999;
  b1_             = -999;
  b2_             = -999;
  q1_             = -999;
  q2_             = -999;
  pt_             ->clear();
  eta_            ->clear();
  phi_            ->clear();
  mass_           ->clear();
  chf_            ->clear();
  nhf_            ->clear();
  phf_            ->clear();
  elf_            ->clear();
  muf_            ->clear();
  jec_            ->clear();
  btag_           ->clear();
  btagIdx_        ->clear();
  etaIdx_         ->clear();
  beta_           ->clear();
  puMva_          ->clear();
  unc_            ->clear();
  ptD_            ->clear();
  ptD_QC_         ->clear();
  vtx3dL_         ->clear();
  vtx3deL_        ->clear();
  vtxPt_          ->clear();
  part_           ->clear();
  puIdL_          ->clear();
  puIdM_          ->clear();
  puIdT_          ->clear();
  idL_            ->clear();
  idM_            ->clear();
  idT_            ->clear();
  axisMinor_      ->clear();
  axisMajor_      ->clear();
  axisMinor_QC_   ->clear();
  axisMajor_QC_   ->clear();
  pull_           ->clear(); 
  pull_QC_        ->clear();
  jetR_           ->clear();
  jetRChg_QC_     ->clear();
  nChg_QC_        ->clear();
  nChg_ptCut_     ->clear();
  nNeutral_ptCut_ ->clear();
  triggerResult_  ->clear();
  softTrackJetPt_ ->clear(); 
  softTrackJetEta_->clear();
  softTrackJetPhi_->clear();
  softTrackJetE_  ->clear();
  //----- MC -------
  npu_ = -999;
  partonDRbb_ = -999;
  partonSt_      ->clear();
  partonId_      ->clear();
  partonMatchIdx_->clear();
  partonMatchDR_ ->clear();
  partonPt_      ->clear();
  partonEta_     ->clear();
  partonPhi_     ->clear();
  partonE_       ->clear();
  genjetPt_      ->clear();
  genjetEta_     ->clear();
  genjetPhi_     ->clear();
  genjetE_       ->clear();
}
//////////////////////////////////////////////////////////////////////////////////////////
VbfHbbFlatTreeProducer::~VbfHbbFlatTreeProducer() 
{
}

DEFINE_FWK_MODULE(VbfHbbFlatTreeProducer);
