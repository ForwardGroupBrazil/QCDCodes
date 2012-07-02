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

#include "KKousour/ZJetAnalysis/plugins/FlatTreeProducer.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "AnalysisDataFormats/CMGTools/interface/Muon.h"
#include "AnalysisDataFormats/CMGTools/interface/Electron.h"
#include "AnalysisDataFormats/CMGTools/interface/Photon.h"
#include "AnalysisDataFormats/CMGTools/interface/DiObject.h"

using namespace std;
using namespace reco;

FlatTreeProducer::FlatTreeProducer(edm::ParameterSet const& cfg) 
{
  srcJets_        = cfg.getParameter<edm::InputTag> ("jets");
  srcMET_         = cfg.getParameter<edm::InputTag> ("met");
  srcElectrons_   = cfg.getParameter<edm::InputTag> ("electrons");
  srcMuons_       = cfg.getParameter<edm::InputTag> ("muons");
  srcDiElectrons_ = cfg.getParameter<edm::InputTag> ("diElectrons");
  srcDiMuons_     = cfg.getParameter<edm::InputTag> ("diMuons");
  srcPhotons_     = cfg.getParameter<edm::InputTag> ("photons");
  srcRho_         = cfg.getParameter<edm::InputTag> ("rho"); 
  minJetPt_       = cfg.getParameter<double>        ("minJetPt");
  maxJetEta_      = cfg.getParameter<double>        ("maxJetEta");   
  minPhotonPt_    = cfg.getParameter<double>        ("minPhotonPt");     
}
//////////////////////////////////////////////////////////////////////////////////////////
void FlatTreeProducer::beginJob() 
{
  jetPt_        = new std::vector<float>(); 
  jetEta_       = new std::vector<float>();
  jetPhi_       = new std::vector<float>();
  jetE_         = new std::vector<float>();
  jetBtag_      = new std::vector<float>();
  jetRMS_       = new std::vector<float>();
  jetBeta_      = new std::vector<float>();
  jetChf_       = new std::vector<float>();
  jetNhf_       = new std::vector<float>();
  jetPhf_       = new std::vector<float>();
  jetMuf_       = new std::vector<float>();
  jetElf_       = new std::vector<float>(); 
  jetId_        = new std::vector<int>();
  jetPu_        = new std::vector<int>();
  photonPt_     = new std::vector<float>();
  photonEta_    = new std::vector<float>();
  photonPhi_    = new std::vector<float>();
  photonE_      = new std::vector<float>();
  elP4_         = new std::vector<LorentzVector>();
  elPt_         = new std::vector<float>();
  elEta_        = new std::vector<float>();
  elPhi_        = new std::vector<float>();
  elE_          = new std::vector<float>(); 
  elIso_        = new std::vector<float>(); 
  elId_         = new std::vector<int>();
  elCh_         = new std::vector<int>();
  muP4_         = new std::vector<LorentzVector>(); 
  muPt_         = new std::vector<float>();
  muEta_        = new std::vector<float>();
  muPhi_        = new std::vector<float>();
  muE_          = new std::vector<float>();
  muIso_        = new std::vector<float>(); 
  muId_         = new std::vector<int>();
  muCh_         = new std::vector<int>();
  
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"       ,&run_         ,"run_/I");
  outTree_->Branch("evtNo"       ,&evt_         ,"evt_/I");
  outTree_->Branch("lumi"        ,&lumi_        ,"lumi_/I");
  outTree_->Branch("nvtx"        ,&nVtx_        ,"nVtx_/I");
  outTree_->Branch("njets"       ,&njets_       ,"njets_/I");
  outTree_->Branch("nmuons"      ,&nmuons_      ,"nmuons_/I");
  outTree_->Branch("nelectrons"  ,&nelectrons_  ,"nelectrons_/I");
  outTree_->Branch("nphotons"    ,&nphotons_    ,"nphotons_/I");
  outTree_->Branch("rho"         ,&rho_         ,"rho_/F");
  outTree_->Branch("met"         ,&met_         ,"met_/F");
  outTree_->Branch("metPhi"      ,&metPhi_      ,"metPhi_/F");
  //---------- jets -----------------------------------------
  outTree_->Branch("jetPt"    ,"vector<float>"     ,&jetPt_);
  outTree_->Branch("jetEta"   ,"vector<float>"     ,&jetEta_);
  outTree_->Branch("jetPhi"   ,"vector<float>"     ,&jetPhi_);
  outTree_->Branch("jetE"     ,"vector<float>"     ,&jetE_);
  outTree_->Branch("jetBtag"  ,"vector<float>"     ,&jetBtag_);
  outTree_->Branch("jetRMS"   ,"vector<float>"     ,&jetRMS_);
  outTree_->Branch("jetBeta"  ,"vector<float>"     ,&jetBeta_);
  outTree_->Branch("jetChf"   ,"vector<float>"     ,&jetChf_);
  outTree_->Branch("jetNhf"   ,"vector<float>"     ,&jetNhf_);
  outTree_->Branch("jetPhf"   ,"vector<float>"     ,&jetPhf_);
  outTree_->Branch("jetMuf"   ,"vector<float>"     ,&jetMuf_);
  outTree_->Branch("jetElf"   ,"vector<float>"     ,&jetElf_);
  outTree_->Branch("jetId"    ,"vector<int>"       ,&jetId_);
  outTree_->Branch("jetPu"    ,"vector<int>"       ,&jetPu_);
  //---------- photons -----------------------------------------
  outTree_->Branch("photonPt" ,"vector<float>"     ,&photonPt_);
  outTree_->Branch("photonEta","vector<float>"     ,&photonEta_);
  outTree_->Branch("photonPhi","vector<float>"     ,&photonPhi_);
  outTree_->Branch("photonE"  ,"vector<float>"     ,&photonE_);
  //---------- electrons -----------------------------------------
  outTree_->Branch("elPt"     ,"vector<float>"     ,&elPt_);
  outTree_->Branch("elEta"    ,"vector<float>"     ,&elEta_);
  outTree_->Branch("elPhi"    ,"vector<float>"     ,&elPhi_);
  outTree_->Branch("elE"      ,"vector<float>"     ,&elE_);
  outTree_->Branch("elIso"    ,"vector<float>"     ,&elIso_);
  outTree_->Branch("elId"     ,"vector<int>"       ,&elId_);
  outTree_->Branch("elCh"     ,"vector<int>"       ,&elCh_);
  //---------- muons -----------------------------------------
  outTree_->Branch("muPt"     ,"vector<float>"     ,&muPt_);
  outTree_->Branch("muEta"    ,"vector<float>"     ,&muEta_);
  outTree_->Branch("muPhi"    ,"vector<float>"     ,&muPhi_);
  outTree_->Branch("muE"      ,"vector<float>"     ,&muE_);
  outTree_->Branch("muIso"    ,"vector<float>"     ,&muIso_);
  outTree_->Branch("muId"     ,"vector<int>"       ,&muId_); 
  outTree_->Branch("muCh"     ,"vector<int>"       ,&muCh_);
  //---------- di-electrons -----------------------------------------
  outTree_->Branch("eePt"     ,&eePt_              ,"eePt_/F");
  outTree_->Branch("eeEta"    ,&eeEta_             ,"eeEta_/F");
  outTree_->Branch("eePhi"    ,&eePhi_             ,"eePhi_/F");
  outTree_->Branch("eeE"      ,&eeE_               ,"eeE_/F");
  outTree_->Branch("eeM"      ,&eeM_               ,"eeM_/F");
  //---------- di-muons -----------------------------------------
  outTree_->Branch("mmPt"     ,&mmPt_              ,"mmPt_/F");
  outTree_->Branch("mmEta"    ,&mmEta_             ,"mmEta_/F");
  outTree_->Branch("mmPhi"    ,&mmPhi_             ,"mmPhi_/F");
  outTree_->Branch("mmE"      ,&mmE_               ,"mmE_/F");
  outTree_->Branch("mmM"      ,&mmM_               ,"mmM_/F");
  //---------- di-electrons+jet -----------------------------------------
  outTree_->Branch("eejPt"    ,&eejPt_             ,"eejPt_/F");
  outTree_->Branch("eejY"     ,&eejY_              ,"eejY_/F");
  outTree_->Branch("eejM"     ,&eejM_              ,"eejM_/F");
  //---------- di-muons+jet -----------------------------------------
  outTree_->Branch("mmjPt"    ,&mmjPt_             ,"mmjPt_/F");
  outTree_->Branch("mmjY"     ,&mmjY_              ,"mmjY_/F");
  outTree_->Branch("mmjM"     ,&mmjM_              ,"mmjM_/F");
  //---------- photon+jet -----------------------------------------
  outTree_->Branch("gjPt"     ,&gjPt_              ,"gjPt_/F");
  outTree_->Branch("gjY"      ,&gjY_               ,"gjY_/F");
  outTree_->Branch("gjM"      ,&gjM_               ,"gjM_/F");
}
//////////////////////////////////////////////////////////////////////////////////////////
void FlatTreeProducer::endJob() 
{
}
//////////////////////////////////////////////////////////////////////////////////////////
void FlatTreeProducer::initialize()
{
  run_          = -999;
  evt_          = -999;
  lumi_         = -999;
  nVtx_         = -999;
  rho_          = -999;
  met_          = -999;
  metPhi_       = -999;
  njets_        = -999;
  nelectrons_   = -999;
  nmuons_       = -999; 
  nphotons_     = -999;
  eePt_         = -999;
  eeEta_        = -999;
  eePhi_        = -999;
  eeE_          = -999;
  eeM_          = -999;
  mmPt_         = -999;
  mmEta_        = -999;
  mmPhi_        = -999;
  mmE_          = -999;
  mmM_          = -999;
  eejPt_        = -999;
  eejY_         = -999;
  eejM_         = -999;
  mmjPt_        = -999;
  mmjY_         = -999;
  mmjM_         = -999; 
  gjPt_         = -999;
  gjY_          = -999;
  gjM_          = -999;
  elP4_->clear();
  elPt_->clear();
  elEta_->clear();
  elPhi_->clear();
  elE_->clear();
  elIso_->clear();
  elId_->clear();
  elCh_->clear();
  muP4_->clear();
  muPt_->clear();
  muEta_->clear();
  muPhi_->clear();
  muE_->clear();
  muIso_->clear();
  muId_->clear();
  muCh_->clear();
  photonPt_->clear();
  photonEta_->clear();
  photonPhi_->clear();
  photonE_->clear();
  jetPt_->clear();
  jetEta_->clear();
  jetPhi_->clear();
  jetE_->clear();
  jetBtag_->clear(); 
  jetRMS_->clear();
  jetBeta_->clear();
  jetChf_->clear();
  jetNhf_->clear();
  jetPhf_->clear();
  jetMuf_->clear();
  jetElf_->clear();
  jetId_->clear();
  jetPu_->clear();
}
//////////////////////////////////////////////////////////////////////////////////////////

void FlatTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  initialize();
 
  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);
  rho_ = *rho;
  //---- reco vertices block --------------------------------------------
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("offlinePrimaryVertices",recVtxs);
  nVtx_ = 0;
  for(VertexCollection::const_iterator i_vtx = recVtxs->begin(); i_vtx != recVtxs->end(); ++i_vtx) {  
    if (!i_vtx->isFake() && (fabs(i_vtx->z()) < 24) && (i_vtx->ndof() >= 4)) {
      nVtx_++;
    }
  }
  //---- met block ---------------------------------------------
  edm::Handle<edm::View<MET> >  met;
  iEvent.getByLabel(srcMET_,met);
  met_ = (*met)[0].et();
  metPhi_ = (*met)[0].phi();
  //---- photon block --------------------------------------------  
  edm::Handle<edm::View<cmg::Photon> > photons;
  iEvent.getByLabel(srcPhotons_,photons);
  edm::View<cmg::Photon> cmg_photons = *photons;
  nphotons_ = 0;
  for(edm::View<cmg::Photon>::const_iterator iph = cmg_photons.begin();iph != cmg_photons.end(); ++iph) { 
    if (iph->pt() < minPhotonPt_) continue;
    photonPt_->push_back(iph->pt());
    photonEta_->push_back(iph->eta());
    photonPhi_->push_back(iph->phi());	
    photonE_->push_back(iph->energy());	 
    nphotons_++;     
  }// electron loop
  //---- muons block --------------------------------------------  
  nmuons_ = 0;
  edm::Handle<edm::View<cmg::Muon> > muons;
  iEvent.getByLabel(srcMuons_,muons);
  edm::View<cmg::Muon> cmg_muons = *muons;
  for(edm::View<cmg::Muon>::const_iterator imuon = cmg_muons.begin();imuon != cmg_muons.end(); ++imuon) { 
    muPt_ ->push_back(imuon->pt());
    muEta_->push_back(imuon->eta());
    muPhi_->push_back(imuon->phi());	
    muE_  ->push_back(imuon->energy());	
    muP4_ ->push_back(imuon->p4());
    int id(0);
    if (
      imuon->isGlobalMuon() && 
      fabs(imuon->normalizedChi2()) < 10 &&
      imuon->numberOfValidMuonHits() > 0 &&
      imuon->numberOfMatches() >1 && 
      fabs(imuon->dxy()) < 0.2 &&
      imuon->numberOfValidPixelHits() > 0 && 
      imuon->numberOfValidTrackerHits() > 10
    ) {
      id = 1;
    }
    muId_->push_back(id);
    muIso_->push_back(-999);
    muCh_->push_back(imuon->charge());      
    nmuons_++;
  }// muon loop 
  //---- electrons block --------------------------------------------  
  nelectrons_ = 0;
  edm::Handle<edm::View<cmg::Electron> > electrons;
  iEvent.getByLabel(srcElectrons_,electrons);
  edm::View<cmg::Electron> cmg_electrons = *electrons;
  for(edm::View<cmg::Electron>::const_iterator iel = cmg_electrons.begin();iel != cmg_electrons.end(); ++iel) { 
    elPt_ ->push_back(iel->pt());
    elEta_->push_back(iel->eta());
    elPhi_->push_back(iel->phi());	
    elE_  ->push_back(iel->energy());
    elP4_ ->push_back(iel->p4());
    float sigmaIetaIeta                  = iel->sigmaIetaIeta();
    float hadronicOverEm                 = iel->hadronicOverEm();
    float deltaPhiSuperClusterTrackAtVtx = iel->deltaPhiSuperClusterTrackAtVtx();
    float deltaEtaSuperClusterTrackAtVtx = iel->deltaEtaSuperClusterTrackAtVtx();
    int id(0);
    float etaSC = fabs(iel->eta());// needs to be FIXED: iel->sourcePtr()->superCluster()->eta()
    if (etaSC < 1.4442) {
      if (sigmaIetaIeta < 0.01 && deltaPhiSuperClusterTrackAtVtx < 0.8 && deltaEtaSuperClusterTrackAtVtx < 0.007 && hadronicOverEm < 0.15) 
      id = 1;
    }// if EB
    if (etaSC > 1.5660) {
      if (sigmaIetaIeta < 0.03 && deltaPhiSuperClusterTrackAtVtx < 0.7 && deltaEtaSuperClusterTrackAtVtx < 0.009 && hadronicOverEm < 0.15) 
      id = 1;
    }// if EE	
    elId_->push_back(id);
    elIso_->push_back(-999);
    elCh_->push_back(iel->charge());
    nelectrons_++;	      
  }// electron loop
  //---- di-electrons block --------------------------------------------  
  edm::Handle<edm::View<cmg::DiObject<cmg::Electron, cmg::Electron> > > dielectrons;
  iEvent.getByLabel(srcDiElectrons_,dielectrons);
  edm::View<cmg::DiObject<cmg::Electron, cmg::Electron> > cmg_dielectrons = *dielectrons;
  if (cmg_dielectrons.size() > 0) {
    eePt_  = cmg_dielectrons[0].pt();
    eeEta_ = cmg_dielectrons[0].eta();
    eePhi_ = cmg_dielectrons[0].phi();
    eeE_   = cmg_dielectrons[0].energy();
    eeM_   = cmg_dielectrons[0].mass();
  }
  //---- di-muon block --------------------------------------------  
  edm::Handle<edm::View<cmg::DiObject<cmg::Muon, cmg::Muon> > > dimuons;
  iEvent.getByLabel(srcDiMuons_,dimuons);
  edm::View<cmg::DiObject<cmg::Muon, cmg::Muon> > cmg_dimuons = *dimuons;
  if (cmg_dimuons.size() > 0) {
    mmPt_  = cmg_dimuons[0].pt();
    mmEta_ = cmg_dimuons[0].eta();
    mmPhi_ = cmg_dimuons[0].phi();
    mmE_   = cmg_dimuons[0].energy();
    mmM_   = cmg_dimuons[0].mass();
  }
  //---- jets block --------------------------------------------  
  njets_ = 0;
  edm::Handle<edm::View<cmg::PFJet> > jets;
  iEvent.getByLabel(srcJets_,jets);
  edm::View<cmg::PFJet> cmg_jets = *jets;
  for(edm::View<cmg::PFJet>::const_iterator ijet = cmg_jets.begin();ijet != cmg_jets.end(); ++ijet) { 
    //---- remove jets beyond acceptance --------
    if (ijet->pt() < minJetPt_ || fabs(ijet->eta()) > maxJetEta_) continue;
    //---- cross clean with electrons and muons --
    bool matched(false);
    for(int iel=0;iel<nelectrons_;iel++) {
      double dR = deltaR(ijet->p4(),elP4_->at(iel));
      if (dR < 0.25) {
        matched = true;
        continue;
      }
    }
    for(int imu=0;imu<nmuons_;imu++) {
      double dR = deltaR(ijet->p4(),muP4_->at(imu));
      if (dR < 0.25) {
        matched = true;
        continue;
      }
    }
    if (matched) continue;

    jetPt_->push_back(ijet->pt());
    jetEta_->push_back(ijet->eta());
    jetPhi_->push_back(ijet->phi());	
    jetE_->push_back(ijet->energy());	
    jetBtag_->push_back(ijet->bDiscriminator("combinedSecondaryVertexBJetTags"));
    jetRMS_->push_back(ijet->rms());
    jetBeta_->push_back(ijet->beta());
    float chf = ijet->component(1).fraction();
    float chm = ijet->component(1).number();
    float nhf = ijet->component(5).fraction();
    float phf = ijet->component(4).fraction();
    float muf = ijet->component(3).fraction();
    float elf = ijet->component(2).fraction();
    int npr = ijet->nConstituents();
    int id(0);
    if (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(ijet->eta())<=2.4 && elf<0.99 && chf>0 && chm>0) || fabs(ijet->eta())>2.4)) {
      id = 1;
    }
    if (npr>1 && phf<0.9 && nhf<0.9 && ((fabs(ijet->eta())<=2.4 && elf<0.9 && muf<0.9 && chf>0 && chm>0) || fabs(ijet->eta())>2.4)) {
      id = 2;
    }
    jetId_->push_back(id);
    // ---placeholder for pu id
    //jetPu_->push_back(ijet->puId("cut-based"));
    jetPu_->push_back(ijet->beta() > 0.1);
    jetChf_->push_back(chf);
    jetNhf_->push_back(nhf);
    jetPhf_->push_back(phf);
    jetMuf_->push_back(muf);
    jetElf_->push_back(elf);	      
    njets_++;
  }// jet loop
  run_    = iEvent.id().run();
  evt_    = iEvent.id().event();
  lumi_   = iEvent.id().luminosityBlock();
  if (nVtx_ > 0 && (nelectrons_ > 1 || nmuons_ > 1 || nphotons_ > 0) && njets_ > 0) {
    if (nelectrons_ > 1) {
      LorentzVector eejP4(cmg_dielectrons[0].p4() + cmg_jets[0].p4()); 
      eejPt_ = eejP4.Pt();
      eejY_  = eejP4.Rapidity();
      eejM_  = eejP4.M();
    }
    if (nmuons_ > 1) {
      LorentzVector mmjP4(cmg_dimuons[0].p4() + cmg_jets[0].p4());
      mmjPt_ = mmjP4.Pt();
      mmjY_  = mmjP4.Rapidity();
      mmjM_  = mmjP4.M();
    }
    if (nphotons_ > 0) { 
      LorentzVector gjP4(cmg_photons[0].p4() + cmg_jets[0].p4());
      gjPt_  = gjP4.Pt();
      gjY_   = gjP4.Rapidity();
      gjM_   = gjP4.M();
    }
    outTree_->Fill();
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
FlatTreeProducer::~FlatTreeProducer() 
{
}

DEFINE_FWK_MODULE(FlatTreeProducer);
