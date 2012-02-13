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

#include "KKousour/MultiJetAnalysis/plugins/PatVBFTree.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

using namespace std;
using namespace reco;

PatVBFTree::PatVBFTree(edm::ParameterSet const& cfg) 
{
  srcJets_ = cfg.getParameter<edm::InputTag>        ("jets");
  srcMET_  = cfg.getParameter<edm::InputTag>        ("met");
  srcRho_  = cfg.getParameter<edm::InputTag>        ("rho");
  srcBtag_ = cfg.getParameter<std::string>          ("btagger");
  srcPU_   = cfg.getUntrackedParameter<std::string> ("pu","");
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatVBFTree::beginJob() 
{
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"       ,&run_         ,"run_/I");
  outTree_->Branch("evtNo"       ,&evt_         ,"evt_/I");
  outTree_->Branch("lumi"        ,&lumi_        ,"lumi_/I");
  outTree_->Branch("nvtx"        ,&nVtx_        ,"nVtx_/I");
  outTree_->Branch("btagIdx"     ,&btagIdx_     ,"btagIdx_[5]/I");
  outTree_->Branch("pvx"         ,&pvx_         ,"pvx_/F");
  outTree_->Branch("pvy"         ,&pvy_         ,"pvy_/F");
  outTree_->Branch("pvz"         ,&pvz_         ,"pvz_/F");
  outTree_->Branch("rho"         ,&rho_         ,"rho_/F");
  outTree_->Branch("ht"          ,&ht_          ,"ht_/F");
  outTree_->Branch("htAll"       ,&htAll_       ,"htAll_/F");
  outTree_->Branch("pvz"         ,&pvz_         ,"pvz_/F");
  outTree_->Branch("met"         ,&met_         ,"met_/F");
  outTree_->Branch("metSig"      ,&metSig_      ,"metSig_/F");
  outTree_->Branch("mqq"         ,&mqq_         ,"mqq_/F");
  outTree_->Branch("mbb"         ,&mbb_         ,"mbb_/F");
  outTree_->Branch("dEtaqq"      ,&dEtaqq_      ,"dEtaqq_/F");
  outTree_->Branch("jetPt"       ,&pt_          ,"pt_[5]/F");
  outTree_->Branch("jetBtag"     ,&btag_        ,"btag_[5]/F");
  outTree_->Branch("jetJec"      ,&jec_         ,"jec_[5]/F");
  outTree_->Branch("jetUnc"      ,&unc_         ,"unc_[5]/F");
  outTree_->Branch("jetBeta"     ,&beta_        ,"beta_[5]/F");
  outTree_->Branch("jetEta"      ,&eta_         ,"eta_[5]/F");
  outTree_->Branch("jetPhi"      ,&phi_         ,"phi_[5]/F");
  outTree_->Branch("jetMass"     ,&mass_        ,"mass_[5]/F");
  outTree_->Branch("jetChf"      ,&chf_         ,"chf_[5]/F");
  outTree_->Branch("jetNhf"      ,&nhf_         ,"nhf_[5]/F");
  outTree_->Branch("jetPhf"      ,&phf_         ,"phf_[5]/F");
  outTree_->Branch("jetMuf"      ,&muf_         ,"muf_[5]/F");
  outTree_->Branch("jetElf"      ,&elf_         ,"elf_[5]/F");
  outTree_->Branch("jetAxis"     ,&axis_        ,"axis_[2][5]/F");
  outTree_->Branch("jetPtD"      ,&ptD_         ,"ptD_[5]/F");
  outTree_->Branch("jetPtMax"    ,&ptMax_       ,"ptMax_[5]/F");
  outTree_->Branch("jetTana"     ,&tana_        ,"tana_[5]/F");
  outTree_->Branch("jetTtheta"   ,&ttheta_      ,"ttheta_[5]/F");
  //------------------- MC ---------------------------------
  outTree_->Branch("inpu"        ,&inpu_        ,"inpu_/F");
  outTree_->Branch("otpu"        ,&otpu_        ,"otpu_/F");
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatVBFTree::endJob() 
{

}
//////////////////////////////////////////////////////////////////////////////////////////
void PatVBFTree::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(srcJets_,jets);
  edm::View<pat::Jet> pat_jets = *jets;

  edm::Handle<edm::View<MET> >  met;
  iEvent.getByLabel(srcMET_,met);

  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("goodOfflinePrimaryVertices",recVtxs);

  edm::Handle<GenEventInfoProduct> hEventInfo;
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;

  initialize();

  //----- at least one good vertex -----------
  bool cut_vtx = (recVtxs->size() > 0);
  //----- at least 4 jets --------------------
  bool cut_njets = (pat_jets.size() > 3);

  if (cut_vtx && cut_njets) {
    bool cutID(true);
    int N(0);
    int Nmax = TMath::Min(int(pat_jets.size()),5);
    float htAll(0.0);
    for(edm::View<pat::Jet>::const_iterator ijet = pat_jets.begin();ijet != pat_jets.end(); ++ijet) { 
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
        chf_[N]     = chf;
        nhf_[N]     = nhf;
        phf_[N]     = phf;
        elf_[N]     = elf;
        muf_[N]     = muf;
        jec_[N]     = 1./ijet->jecFactor(0);
        pt_[N]      = ijet->pt();
        phi_[N]     = ijet->phi();
        eta_[N]     = ijet->eta();
        mass_[N]    = ijet->mass();
        btag_[N]    = ijet->bDiscriminator(srcBtag_);
        beta_[N]    = ijet->userFloat("beta");
        unc_[N]     = ijet->userFloat("jecUnc");
        ptD_[N]     = ijet->userFloat("ptD");
        ptMax_[N]   = ijet->userFloat("ptMax");
        axis_[0][N] = ijet->userFloat("axis1");
        axis_[1][N] = ijet->userFloat("axis2");
        tana_[N]    = ijet->userFloat("tana");
        ttheta_[N]  = ijet->userFloat("ttheta");
      }
      N++;
    }// jet loop
    //----- order jets according to the btag value -----
    float tmp_btag[4],ht(0.0);
    for(int i=0;i<4;i++) {
      btagIdx_[i] = i;
      tmp_btag[i] = btag_[i];
      ht += pt_[i]; 
    }
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
    // ------- MC -----------------------------------------
    int ootpuEarly(0),ootpuLate(0),intpu(0);
    if (!iEvent.isRealData() && srcPU_ != "") {
      iEvent.getByLabel(srcPU_,PupInfo);
      std::vector<PileupSummaryInfo>::const_iterator PUI;
      for(PUI = PupInfo->begin(); PUI != PupInfo->end(); ++PUI) {
        if (PUI->getBunchCrossing() < 0)
          ootpuEarly += PUI->getPU_NumInteractions();
        else if (PUI->getBunchCrossing() > 0)
          ootpuLate += PUI->getPU_NumInteractions();
        else
          intpu += PUI->getPU_NumInteractions();
      }
    }
    if (cutID) {
      // ---- MC ----------
      inpu_   = intpu;
      otpu_   = ootpuEarly + ootpuLate;
      // ------------------
      ht_     = ht;
      htAll_  = htAll;
      rho_    = *rho;
      pvx_    = (*recVtxs)[0].x();
      pvy_    = (*recVtxs)[0].y();
      pvz_    = (*recVtxs)[0].z();
      nVtx_   = recVtxs->size();
      met_    = (*met)[0].et();
      metSig_ = (*met)[0].et()/(*met)[0].sumEt();
      run_    = iEvent.id().run();
      evt_    = iEvent.id().event();
      lumi_   = iEvent.id().luminosityBlock();
      mqq_    = (pat_jets[btagIdx_[2]].p4()+pat_jets[btagIdx_[3]].p4()).mass();
      mbb_    = (pat_jets[btagIdx_[0]].p4()+pat_jets[btagIdx_[1]].p4()).mass();
      dEtaqq_ = fabs(pat_jets[btagIdx_[2]].eta()-pat_jets[btagIdx_[3]].eta());
      outTree_->Fill();
    }
  }// if vtx and jet multi
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatVBFTree::initialize()
{
  run_    = -999;
  evt_    = -999;
  lumi_   = -999;
  nVtx_   = -999;
  rho_    = -999;
  met_    = -999;
  metSig_ = -999;
  ht_     = -999;
  htAll_  = -999;
  pvx_    = -999;
  pvy_    = -999;
  pvz_    = -999;
  mqq_    = -999;
  mbb_    = -999;
  dEtaqq_ = -999;
  for(int i=0;i<5;i++) {
    pt_[i]      = -999;
    eta_[i]     = -999;
    phi_[i]     = -999;
    mass_[i]    = -999;
    chf_[i]     = -999;
    nhf_[i]     = -999;
    phf_[i]     = -999;
    elf_[i]     = -999;
    muf_[i]     = -999;
    jec_[i]     = -999;
    btag_[i]    = -999;
    beta_[i]    = -999;
    unc_[i]     = -999;
    ptD_[i]     = -999;
    ptMax_[i]   = -999;
    axis_[0][i] = -999;
    axis_[1][i] = -999;
    tana_[i]    = -999;
    ttheta_[i]  = -999;
  }  
}
//////////////////////////////////////////////////////////////////////////////////////////
PatVBFTree::~PatVBFTree() 
{
}

DEFINE_FWK_MODULE(PatVBFTree);
