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
#include "TLorentzVector.h"

#include "KKousour/MultiJetAnalysis/plugins/PatVBFTree.h"
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

using namespace std;
using namespace reco;

PatVBFTree::PatVBFTree(edm::ParameterSet const& cfg) 
{
  srcJets_            = cfg.getParameter<edm::InputTag>             ("jets");
  srcMET_             = cfg.getParameter<edm::InputTag>             ("met");
  srcRho_             = cfg.getParameter<edm::InputTag>             ("rho");
  srcRhoQGL_          = cfg.getParameter<edm::InputTag>             ("rhoQGL"); 
  srcBtag_            = cfg.getParameter<std::string>               ("btagger");
  dEtaMin_            = cfg.getParameter<double>                    ("dEtaMin");
  ptMin_              = cfg.getParameter<vector<double> >           ("ptMin");
  srcPU_              = cfg.getUntrackedParameter<std::string>      ("pu","");
  srcGenJets_         = cfg.getUntrackedParameter<edm::InputTag>    ("genjets",edm::InputTag(""));

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
void PatVBFTree::beginJob() 
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
  outTree_->Branch("btagIdx"              ,&btagIdx_           ,"btagIdx_[5]/I");
  outTree_->Branch("etaIdx"               ,&etaIdx_            ,"etaIdx_[5]/I");
  outTree_->Branch("softHt"               ,&softHt_            ,"softHt_/F");
  outTree_->Branch("pvx"                  ,&pvx_               ,"pvx_/F");
  outTree_->Branch("pvy"                  ,&pvy_               ,"pvy_/F");
  outTree_->Branch("pvz"                  ,&pvz_               ,"pvz_/F");
  outTree_->Branch("rho"                  ,&rho_               ,"rho_/F");
  outTree_->Branch("rhoQGL"               ,&rhoQGL_            ,"rhoQGL_/F");
  outTree_->Branch("ht"                   ,&ht_                ,"ht_/F");
  outTree_->Branch("htAll"                ,&htAll_             ,"htAll_/F");
  outTree_->Branch("pvz"                  ,&pvz_               ,"pvz_/F");
  outTree_->Branch("met"                  ,&met_               ,"met_/F");
  outTree_->Branch("metPhi"               ,&metPhi_            ,"metPhi_/F");
  outTree_->Branch("metSig"               ,&metSig_            ,"metSig_/F");
  outTree_->Branch("mqq"                  ,&mqq_               ,"mqq_/F");
  outTree_->Branch("mbb"                  ,&mbb_               ,"mbb_/F");
  outTree_->Branch("dEtaqq"               ,&dEtaqq_            ,"dEtaqq_/F");
  outTree_->Branch("mqqEta"               ,&mqqEta_            ,"mqqEta_/F");
  outTree_->Branch("mbbEta"               ,&mbbEta_            ,"mbbEta_/F");
  outTree_->Branch("dEtaqqEta"            ,&dEtaqqEta_         ,"dEtaqqEta_/F");
  outTree_->Branch("dEtabb"               ,&dEtabb_            ,"dEtabb_/F");
  outTree_->Branch("ptqq"                 ,&ptqq_              ,"ptqq_/F");
  outTree_->Branch("ptbb"                 ,&ptbb_              ,"ptbb_/F");
  outTree_->Branch("dPhiqq"               ,&dPhiqq_            ,"dPhiqq_/F");
  outTree_->Branch("dPhibb"               ,&dPhibb_            ,"dPhibb_/F");
  outTree_->Branch("etaBoostqq"           ,&etaBoostqq_        ,"etaBoostqq_/F");
  outTree_->Branch("etaBoostbb"           ,&etaBoostbb_        ,"etaBoostbb_/F");
  outTree_->Branch("jetPt"                ,&pt_                ,"pt_[5]/F");
  outTree_->Branch("jetBtag"              ,&btag_              ,"btag_[5]/F");
  outTree_->Branch("jetJec"               ,&jec_               ,"jec_[5]/F");
  outTree_->Branch("jetUnc"               ,&unc_               ,"unc_[5]/F");
  outTree_->Branch("jetBeta"              ,&beta_              ,"beta_[5]/F");
  outTree_->Branch("jetQGL"               ,&qgl_               ,"qgl_[5]/F");
  outTree_->Branch("jetEta"               ,&eta_               ,"eta_[5]/F");
  outTree_->Branch("jetPhi"               ,&phi_               ,"phi_[5]/F");
  outTree_->Branch("jetMass"              ,&mass_              ,"mass_[5]/F");
  outTree_->Branch("jetChf"               ,&chf_               ,"chf_[5]/F");
  outTree_->Branch("jetNhf"               ,&nhf_               ,"nhf_[5]/F");
  outTree_->Branch("jetPhf"               ,&phf_               ,"phf_[5]/F");
  outTree_->Branch("jetMuf"               ,&muf_               ,"muf_[5]/F");
  outTree_->Branch("jetElf"               ,&elf_               ,"elf_[5]/F");
  outTree_->Branch("jetPtD"               ,&ptD_               ,"ptD_[5]/F");
 //------------------------------------------------------------------- 
  outTree_->Branch("jetAxis"              ,&axis_              ,"axis_[2][5]/F");
  outTree_->Branch("jetAxis_QC"           ,&axis_QC_           ,"axis_QC_[2][5]/F");
  outTree_->Branch("jetPull"              ,&pull_              ,"pull_[5]/F");
  outTree_->Branch("jetPull_QC"           ,&pull_QC_           ,"pull_QC_[5]/F");
  outTree_->Branch("jetR"                 ,&jetR_              ,"jetR_[5]/F");
  outTree_->Branch("jetRChg_QC"           ,&jetRChg_QC_        ,"jetRChg_QC_[5]/F");
  outTree_->Branch("jetnChg_QC"           ,&nChg_QC_           ,"nChg_QC_[5]/I");
  outTree_->Branch("jetnChg_ptCut"        ,&nChg_ptCut_        ,"nChg_ptCut_[5]/I");
  outTree_->Branch("jetnNeutral_ptCut"    ,&nNeutral_ptCut_    ,"nNeutral_ptCut_[5]/I");
  outTree_->Branch("jetVtxMass"           ,&vtxMass_           ,"vtxMass_[5]/F"); 
  outTree_->Branch("jetVtx3dL"            ,&vtx3dL_            ,"vtx3dL_[5]/F");
  outTree_->Branch("jetVtx3deL"           ,&vtx3deL_           ,"vtx3deL_[5]/F");
  outTree_->Branch("jetSumTrkPt"          ,&sumTrkPt_          ,"sumTrkPt_[5]/F");
  outTree_->Branch("jetSumTrkP"           ,&sumTrkP_           ,"sumTrkP_[5]/F");
  outTree_->Branch("jetSumTrkPtV"         ,&sumTrkPtV_         ,"sumTrkPtV_[5]/F");
  outTree_->Branch("jetLeadTrkPt"         ,&leadTrkPt_         ,"leadTrkPt_[5]/F");
  outTree_->Branch("jetVtxPt"             ,&vtxPt_             ,"vtxPt_[5]/F");
  outTree_->Branch("jetVtxNTrks"          ,&vtxNTrks_          ,"vtxNTrks_[5]/I"); 
  outTree_->Branch("jetPart"              ,&part_              ,"part_[5]/I"); 
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
  partonId_  = new std::vector<int>;
  partonSt_  = new std::vector<int>;
  partonPt_  = new std::vector<float>;
  partonEta_ = new std::vector<float>;
  partonPhi_ = new std::vector<float>;
  partonE_   = new std::vector<float>;
  genjetPt_  = new std::vector<float>;
  genjetEta_ = new std::vector<float>;
  genjetPhi_ = new std::vector<float>;
  genjetE_   = new std::vector<float>;
  outTree_->Branch("partonId"       ,"vector<int>"   ,&partonId_);
  outTree_->Branch("partonSt"       ,"vector<int>"   ,&partonSt_);
  outTree_->Branch("partonPt"       ,"vector<float>" ,&partonPt_);
  outTree_->Branch("partonEta"      ,"vector<float>" ,&partonEta_);
  outTree_->Branch("partonPhi"      ,"vector<float>" ,&partonPhi_);
  outTree_->Branch("partonE"        ,"vector<float>" ,&partonE_);
  outTree_->Branch("genjetPt"       ,"vector<float>" ,&genjetPt_);
  outTree_->Branch("genjetEta"      ,"vector<float>" ,&genjetEta_);
  outTree_->Branch("genjetPhi"      ,"vector<float>" ,&genjetPhi_);
  outTree_->Branch("genjetE"        ,"vector<float>" ,&genjetE_);

  cout<<"Begin job finished"<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatVBFTree::endJob() 
{
  delete partonSt_;
  delete partonId_;
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
  for(unsigned i=0;i<vtriggerSelector_.size();i++) {
    delete vtriggerSelector_[i];
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatVBFTree::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByLabel(srcJets_,jets); 
  pat::JetCollection pat_jets = *jets;

  edm::Handle<edm::ValueMap<float> > qglMap;
  iEvent.getByLabel("qglAK5PF",qglMap);

  edm::Handle<reco::TrackJetCollection> softTrackJets;
  iEvent.getByLabel("ak5SoftTrackJets",softTrackJets);

  edm::Handle<edm::View<MET> >  met;
  iEvent.getByLabel(srcMET_,met);

  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);

  edm::Handle<double> rhoQGL;
  iEvent.getByLabel(srcRhoQGL_,rhoQGL);

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("goodOfflinePrimaryVertices",recVtxs);

  edm::Handle<GenEventInfoProduct> hEventInfo;
  edm::Handle<GenParticleCollection> genParticles;
  edm::Handle<GenJetCollection> genjets;
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
  //----- at least 4 jets --------------------
  bool cut_njets = (jets->size() > 3);

  if (cut_vtx && cut_njets) {
    //----- soft track jets ----------------------
    nSoftTrackJets_ = int(softTrackJets->size());
    float softHt(0.0);
    for(int it=0;it<nSoftTrackJets_;it++) {
      softHt += (*softTrackJets)[it].pt();
      if (it < 3) {
        softTrackJetPt_ ->push_back((*softTrackJets)[it].pt());
        softTrackJetEta_->push_back((*softTrackJets)[it].eta());
        softTrackJetPhi_->push_back((*softTrackJets)[it].phi());
        softTrackJetE_  ->push_back((*softTrackJets)[it].energy());
      } 
    }
    //----- PF jets ------------------------------
    bool cutID(true);
    int N(0);
    int Nmax = TMath::Min(int(jets->size()),5);
    float htAll(0.0);
    for(pat::JetCollection::const_iterator ijet = jets->begin();ijet != jets->end(); ++ijet) { 
      bool id(false);
      float chf = ijet->chargedHadronEnergyFraction();
      float nhf = ijet->neutralHadronEnergyFraction() + ijet->HFHadronEnergyFraction();
      float phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      float elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());
      float muf = ijet->muonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      int chm = ijet->chargedHadronMultiplicity();
      int npr = ijet->chargedMultiplicity() + ijet->neutralMultiplicity();
      id = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(ijet->eta())<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || fabs(ijet->eta())>2.4));
      htAll += ijet->pt();
      if (N < Nmax) {
        cutID *= id;
        chf_[N]               = chf;
        nhf_[N]               = nhf;
        phf_[N]               = phf;
        elf_[N]               = elf;
        muf_[N]               = muf;
        jec_[N]               = 1./ijet->jecFactor(0);
        pt_[N]                = ijet->pt();
        phi_[N]               = ijet->phi();
        eta_[N]               = ijet->eta();
        mass_[N]              = ijet->mass();
        btag_[N]              = ijet->bDiscriminator(srcBtag_);
        beta_[N]              = ijet->userFloat("beta");
        unc_[N]               = ijet->userFloat("jecUnc");
        ptD_[N]               = ijet->userFloat("ptD");
        vtxMass_[N]           = ijet->userFloat("vtxMass");
        vtx3dL_[N]            = ijet->userFloat("vtx3dL");
        vtx3deL_[N]           = ijet->userFloat("vtx3deL");
        sumTrkPt_[N]          = ijet->userFloat("sumTrkPt");
        sumTrkP_[N]           = ijet->userFloat("sumTrkP");
        sumTrkPtV_[N]         = ijet->userFloat("sumTrkPtV");
        leadTrkPt_[N]         = ijet->userFloat("leadTrkPt");
        vtxPt_[N]             = ijet->userFloat("vtxPt");
        vtxNTrks_[N]          = ijet->userInt("vtxNTracks");
        part_[N]              = ijet->userInt("nConstituents");

        axis_[0][N]           = ijet->userFloat("axis1");
        axis_[1][N]           = ijet->userFloat("axis2");
        axis_QC_[0][N]        = ijet->userFloat("axis1_QC");
        axis_QC_[1][N]        = ijet->userFloat("axis2_QC");
        pull_[N]              = ijet->userFloat("pull");
        pull_QC_[N]           = ijet->userFloat("pull_QC");
        jetR_[N]              = ijet->userFloat("jetR");
        jetRChg_QC_[N]        = ijet->userFloat("jetRchg_QC");
        nChg_ptCut_[N]        = ijet->userInt("nChg_ptCut");
        nNeutral_ptCut_[N]    = ijet->userInt("nNeutral_ptCut");
        nChg_QC_[N]           = ijet->userInt("nChg_QC");
        //----- calculate the qg likelihood ------------
        int ii = ijet - jets->begin();
        edm::RefToBase<pat::Jet> jetRef(edm::Ref<pat::JetCollection>(jets,ii));
        qgl_[N] = (*qglMap)[jetRef];
      }
      N++;
    }// jet loop
    //----- order the 4 leading jets according to the btag value -----
    float tmp_btag[4],ht(0.0);
    for(int i=0;i<4;i++) {
      btagIdx_[i] = i;
      tmp_btag[i] = btag_[i];
      ht += pt_[i]; 
    }
    btagIdx_[4] = 4;
    for(int i=0;i<4;i++) {
      for(int j=i+1;j<4;j++) {
        if (tmp_btag[j] > tmp_btag[i]) {
          int k = btagIdx_[i];
          btagIdx_[i] = btagIdx_[j];
          btagIdx_[j] = k;
          float tmp = tmp_btag[i]; 
          tmp_btag[i] = tmp_btag[j];
          tmp_btag[j] = tmp;
        }
      } 
    }
    //----- order the 4 leading jets according to the eta value -----
    float tmp_eta[4];
    for(int i=0;i<4;i++) {
      etaIdx_[i] = i;
      tmp_eta[i] = eta_[i]; 
    }
    etaIdx_[4] = 4;
    for(int i=0;i<4;i++) {
      for(int j=i+1;j<4;j++) {
        if (tmp_eta[j] > tmp_eta[i]) {
          int k = etaIdx_[i];
          etaIdx_[i] = etaIdx_[j];
          etaIdx_[j] = k;
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
      for(reco::GenJetCollection::const_iterator igen = genjets->begin();igen != genjets->end(); ++igen) {
        genjetPt_ ->push_back(igen->pt());
        genjetEta_->push_back(igen->eta());
        genjetPhi_->push_back(igen->phi());
        genjetE_  ->push_back(igen->energy()); 
      }// genjet loop
      //---------- partons ------------------
      iEvent.getByLabel("genParticles", genParticles);
      for(unsigned ip = 0; ip < genParticles->size(); ++ ip) {
        const GenParticle &p = (*genParticles)[ip];
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
      }// parton loop
    }// if MC    
    ht_     = ht;
    htAll_  = htAll;
    softHt_ = softHt;
    rho_    = *rho;
    rhoQGL_ = *rhoQGL;
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
    mqq_    = (pat_jets[btagIdx_[2]].p4()+pat_jets[btagIdx_[3]].p4()).mass();
    mbb_    = (pat_jets[btagIdx_[0]].p4()+pat_jets[btagIdx_[1]].p4()).mass();
    ptqq_   = (pat_jets[btagIdx_[2]].p4()+pat_jets[btagIdx_[3]].p4()).pt();
    ptbb_   = (pat_jets[btagIdx_[0]].p4()+pat_jets[btagIdx_[1]].p4()).pt();
    dEtabb_ = fabs(pat_jets[btagIdx_[0]].eta()-pat_jets[btagIdx_[1]].eta());
    dEtaqq_ = fabs(pat_jets[btagIdx_[2]].eta()-pat_jets[btagIdx_[3]].eta());
    etaBoostbb_ = 0.5*(pat_jets[btagIdx_[0]].eta()+pat_jets[btagIdx_[1]].eta());
    etaBoostqq_ = 0.5*(pat_jets[btagIdx_[2]].eta()+pat_jets[btagIdx_[3]].eta());
    dPhibb_ = fabs(deltaPhi(pat_jets[btagIdx_[0]].phi(),pat_jets[btagIdx_[1]].phi()));
    dPhiqq_ = fabs(deltaPhi(pat_jets[btagIdx_[2]].phi(),pat_jets[btagIdx_[3]].phi()));
    //--- eta ordered quantities ----------
    mqqEta_    = (pat_jets[etaIdx_[0]].p4()+pat_jets[etaIdx_[3]].p4()).mass();
    mbbEta_    = (pat_jets[etaIdx_[1]].p4()+pat_jets[etaIdx_[2]].p4()).mass();
    dEtaqqEta_ = fabs(pat_jets[etaIdx_[0]].eta()-pat_jets[etaIdx_[3]].eta());
    //--- pt cut --------------------------
    bool cut_PT(true);
    for(int j=0;j<TMath::Min(4,int(ptMin_.size()));j++) {
      cut_PT *= (pt_[j] > ptMin_[j]); 
    }
    if (cutID && cut_PT && (dEtaqq_ > dEtaMin_)) {
      outTree_->Fill();
    }
  }// if vtx and jet multi
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatVBFTree::initialize()
{
  run_            = -999;
  evt_            = -999;
  lumi_           = -999;
  nVtx_           = -999;
  nSoftTrackJets_ = -999;
  rho_            = -999;
  rhoQGL_         = -999;
  met_            = -999;
  metPhi_         = -999;
  metSig_         = -999;
  ht_             = -999;
  htAll_          = -999;
  softHt_         = -999;
  pvx_            = -999;
  pvy_            = -999;
  pvz_            = -999;
  mqq_            = -999;
  mbb_            = -999;
  mqqEta_         = -999;
  mbbEta_         = -999;
  ptqq_           = -999;
  ptbb_           = -999;
  dEtaqq_         = -999;
  dEtaqqEta_      = -999;
  dEtabb_         = -999;
  etaBoostqq_     = -999;
  etaBoostbb_     = -999;
  dPhiqq_         = -999;
  dPhibb_         = -999;
  for(int i=0;i<5;i++) {
    pt_[i]        = -999;
    eta_[i]       = -999;
    phi_[i]       = -999;
    mass_[i]      = -999;
    chf_[i]       = -999;
    nhf_[i]       = -999;
    phf_[i]       = -999;
    elf_[i]       = -999;
    muf_[i]       = -999;
    jec_[i]       = -999;
    qgl_[i]       = -999;
    btag_[i]      = -999;
    btagIdx_[i]   = -999;
    etaIdx_[i]    = -999;
    beta_[i]      = -999;
    unc_[i]       = -999;
    ptD_[i]       = -999;
    vtxMass_[i]   = -999;
    vtx3dL_[i]    = -999;
    vtx3deL_[i]   = -999;
    sumTrkPt_[i]  = -999;
    sumTrkP_[i]   = -999;
    sumTrkPtV_[i] = -999;
    leadTrkPt_[i] = -999;
    vtxPt_[i]     = -999;
    vtxNTrks_[i]  = -999;
    part_[i]      = -999;
    for(int j=0;j<2;j++) {
      axis_[j][i]   = -999;
      axis_QC_[j][i] = -999;
    }
    pull_[i]           = -999; 
    pull_QC_[i]        = -999;
    jetR_[i]           = -999;
    jetRChg_QC_[i]     = -999;
    nChg_QC_[i]        = -999;
    nChg_ptCut_[i]     = -999;
    nNeutral_ptCut_[i] = -999;
  }
  triggerResult_  ->clear();
  softTrackJetPt_ ->clear(); 
  softTrackJetEta_->clear();
  softTrackJetPhi_->clear();
  softTrackJetE_  ->clear();
  //----- MC -------
  npu_ = -999;
  partonSt_ ->clear();
  partonId_ ->clear();
  partonPt_ ->clear();
  partonEta_->clear();
  partonPhi_->clear();
  partonE_  ->clear();
  genjetPt_ ->clear();
  genjetEta_->clear();
  genjetPhi_->clear();
  genjetE_  ->clear();
}
//////////////////////////////////////////////////////////////////////////////////////////
PatVBFTree::~PatVBFTree() 
{
}

DEFINE_FWK_MODULE(PatVBFTree);
