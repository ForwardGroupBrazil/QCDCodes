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

using namespace std;
using namespace reco;

VbfHbbFlatTreeProducer::VbfHbbFlatTreeProducer(edm::ParameterSet const& cfg) 
{
  srcJets_            = cfg.getParameter<edm::InputTag>             ("jets");
  srcMET_             = cfg.getParameter<edm::InputTag>             ("met");
  srcRho_             = cfg.getParameter<edm::InputTag>             ("rho"); 
  srcBtag_            = cfg.getParameter<std::string>               ("btagger");
  shiftJES_           = cfg.getParameter<double>                    ("shiftJES");
  dEtaMin_            = cfg.getParameter<double>                    ("dEtaMin");
  ptMin_              = cfg.getParameter<vector<double> >           ("ptMin");
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
  outTree_->Branch("btagIdx"              ,&btagIdx_           ,"btagIdx_[5]/I");
  outTree_->Branch("etaIdx"               ,&etaIdx_            ,"etaIdx_[5]/I");
  outTree_->Branch("softHt"               ,&softHt_            ,"softHt_/F");
  outTree_->Branch("pvx"                  ,&pvx_               ,"pvx_/F");
  outTree_->Branch("pvy"                  ,&pvy_               ,"pvy_/F");
  outTree_->Branch("pvz"                  ,&pvz_               ,"pvz_/F");
  outTree_->Branch("rho"                  ,&rho_               ,"rho_/F");
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
  outTree_->Branch("jetPuMva"             ,&puMva_             ,"puMva_[5]/F");
  outTree_->Branch("jetJec"               ,&jec_               ,"jec_[5]/F");
  outTree_->Branch("jetUnc"               ,&unc_               ,"unc_[5]/F");
  outTree_->Branch("jetBeta"              ,&beta_              ,"beta_[5]/F");
  outTree_->Branch("jetEta"               ,&eta_               ,"eta_[5]/F");
  outTree_->Branch("jetPhi"               ,&phi_               ,"phi_[5]/F");
  outTree_->Branch("jetMass"              ,&mass_              ,"mass_[5]/F");
  outTree_->Branch("jetChf"               ,&chf_               ,"chf_[5]/F");
  outTree_->Branch("jetNhf"               ,&nhf_               ,"nhf_[5]/F");
  outTree_->Branch("jetPhf"               ,&phf_               ,"phf_[5]/F");
  outTree_->Branch("jetMuf"               ,&muf_               ,"muf_[5]/F");
  outTree_->Branch("jetElf"               ,&elf_               ,"elf_[5]/F");
  outTree_->Branch("jetPtD"               ,&ptD_               ,"ptD_[5]/F");
  outTree_->Branch("jetPtD_QC"            ,&ptD_QC_            ,"ptD_QC_[5]/F");
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
  outTree_->Branch("jetVtx3dL"            ,&vtx3dL_            ,"vtx3dL_[5]/F");
  outTree_->Branch("jetVtx3deL"           ,&vtx3deL_           ,"vtx3deL_[5]/F");
  outTree_->Branch("jetVtxPt"             ,&vtxPt_             ,"vtxPt_[5]/F");
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
void VbfHbbFlatTreeProducer::endJob() 
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
void VbfHbbFlatTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  edm::Handle<edm::View<cmg::PFJet> > jets;
  iEvent.getByLabel(srcJets_,jets);
  edm::View<cmg::PFJet> cmg_jets = *jets;

  edm::Handle<edm::ValueMap<float> > qglMap;
  iEvent.getByLabel("qglAK5PF",qglMap);

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
  //edm::Handle<GenJetCollection> genjets;
  edm::Handle<std::vector<cmg::PhysicsObjectWithPtr<edm::Ptr<reco::GenJet> > > > genjets;
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
    bool cutID(true);
    int N(0);
    int Nmax = TMath::Min(int(jets->size()),5);
    float htAll(0.0);
    for(edm::View<cmg::PFJet>::const_iterator ijet = cmg_jets.begin();ijet != cmg_jets.end(); ++ijet) {  
      bool id(false);

      float chf = ijet->component(1).fraction();
      float nhf = ijet->component(5).fraction();
      float phf = ijet->component(4).fraction();
      float muf = ijet->component(3).fraction();
      float elf = ijet->component(2).fraction();
      int chm = ijet->component(1).number();
      int npr = ijet->nConstituents();

      id = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(ijet->eta())<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || fabs(ijet->eta())>2.4));
      htAll += ijet->pt();
      if (N < Nmax) {
        cutID *= id;
        chf_[N]               = chf;
        nhf_[N]               = nhf;
        phf_[N]               = phf;
        elf_[N]               = elf;
        muf_[N]               = muf;
        unc_[N]               = ijet->uncOnFourVectorScale();
        double aJES = (1.+shiftJES_*ijet->uncOnFourVectorScale());
        jec_[N]               = aJES/ijet->rawFactor();
        pt_[N]                = aJES*ijet->pt();
        phi_[N]               = ijet->phi();
        eta_[N]               = ijet->eta();
        mass_[N]              = aJES*ijet->mass();
        btag_[N]              = ijet->bDiscriminator(srcBtag_.c_str());
        beta_[N]              = ijet->beta();
        puMva_[N]             = ijet->puMva("full");
        ptD_[N]               = ijet->ptd();
        ptD_QC_[N]            = ijet->ptdQC();
        vtx3dL_[N]            = ijet->vtx3dL();
        vtx3deL_[N]           = ijet->vtx3deL();
        vtxPt_[N]             = ijet->vtxPt();
        part_[N]              = npr;

        axis_[0][N]           = ijet->axisMajor();
        axis_[1][N]           = ijet->axisMinor();
        axis_QC_[0][N]        = ijet->axisMajorQC();
        axis_QC_[1][N]        = ijet->axisMinorQC();
        pull_[N]              = ijet->pull();
        pull_QC_[N]           = ijet->pullQC();
        jetR_[N]              = ijet->fmax();
        jetRChg_QC_[N]        = ijet->fmaxCharged();
        nChg_ptCut_[N]        = ijet->nChargedPtCut();
        nNeutral_ptCut_[N]    = ijet->nNeutralPtCut();
        nChg_QC_[N]           = ijet->nChargedQC();
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
    double vaJES[4];
    for(int i=0;i<4;i++) {
      vaJES[i] = 1.+shiftJES_*cmg_jets[i].uncOnFourVectorScale();
    }
    mqq_    = (vaJES[btagIdx_[2]]*cmg_jets[btagIdx_[2]].p4()+vaJES[btagIdx_[3]]*cmg_jets[btagIdx_[3]].p4()).mass();
    mbb_    = (vaJES[btagIdx_[0]]*cmg_jets[btagIdx_[0]].p4()+vaJES[btagIdx_[1]]*cmg_jets[btagIdx_[1]].p4()).mass();
    ptqq_   = (vaJES[btagIdx_[2]]*cmg_jets[btagIdx_[2]].p4()+vaJES[btagIdx_[3]]*cmg_jets[btagIdx_[3]].p4()).pt();
    ptbb_   = (vaJES[btagIdx_[0]]*cmg_jets[btagIdx_[0]].p4()+vaJES[btagIdx_[1]]*cmg_jets[btagIdx_[1]].p4()).pt();
    dEtabb_ = fabs(cmg_jets[btagIdx_[0]].eta()-cmg_jets[btagIdx_[1]].eta());
    dEtaqq_ = fabs(cmg_jets[btagIdx_[2]].eta()-cmg_jets[btagIdx_[3]].eta());
    etaBoostbb_ = 0.5*(cmg_jets[btagIdx_[0]].eta()+cmg_jets[btagIdx_[1]].eta());
    etaBoostqq_ = 0.5*(cmg_jets[btagIdx_[2]].eta()+cmg_jets[btagIdx_[3]].eta());
    dPhibb_ = fabs(deltaPhi(cmg_jets[btagIdx_[0]].phi(),cmg_jets[btagIdx_[1]].phi()));
    dPhiqq_ = fabs(deltaPhi(cmg_jets[btagIdx_[2]].phi(),cmg_jets[btagIdx_[3]].phi()));
    //--- eta ordered quantities ----------
    mqqEta_    = (vaJES[etaIdx_[0]]*cmg_jets[etaIdx_[0]].p4()+vaJES[etaIdx_[3]]*cmg_jets[etaIdx_[3]].p4()).mass();
    mbbEta_    = (vaJES[etaIdx_[1]]*cmg_jets[etaIdx_[1]].p4()+vaJES[etaIdx_[2]]*cmg_jets[etaIdx_[2]].p4()).mass();
    dEtaqqEta_ = fabs(cmg_jets[etaIdx_[0]].eta()-cmg_jets[etaIdx_[3]].eta());
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
void VbfHbbFlatTreeProducer::initialize()
{
  run_            = -999;
  evt_            = -999;
  lumi_           = -999;
  nVtx_           = -999;
  nSoftTrackJets_ = -999;
  rho_            = -999;
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
    btag_[i]      = -999;
    btagIdx_[i]   = -999;
    etaIdx_[i]    = -999;
    beta_[i]      = -999;
    puMva_[i]     = -999;
    unc_[i]       = -999;
    ptD_[i]       = -999;
    ptD_QC_[i]    = -999;
    vtx3dL_[i]    = -999;
    vtx3deL_[i]   = -999;
    vtxPt_[i]     = -999;
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
VbfHbbFlatTreeProducer::~VbfHbbFlatTreeProducer() 
{
}

DEFINE_FWK_MODULE(VbfHbbFlatTreeProducer);
