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

#include "KKousour/MultiJetAnalysis/plugins/LowerMultiplicityTree.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

using namespace std;
using namespace reco;

LowerMultiplicityTree::LowerMultiplicityTree(edm::ParameterSet const& cfg) 
{
  srcJets_    = cfg.getParameter<edm::InputTag>          ("jets");
  srcMET_     = cfg.getParameter<edm::InputTag>          ("met");
  srcRho_     = cfg.getParameter<edm::InputTag>          ("rho");
  srcGenJets_ = cfg.getUntrackedParameter<edm::InputTag> ("genjets",edm::InputTag("ak5GenJets"));
  srcPU_      = cfg.getUntrackedParameter<string>        ("pu","");
  isMC_       = cfg.getUntrackedParameter<bool>          ("isMC",false);
  etaMax_     = cfg.getParameter<double>                 ("etaMAX");
  ptMin_      = cfg.getParameter<double>                 ("ptMIN");
  betaMin_    = cfg.getParameter<double>                 ("betaMIN");
}
//////////////////////////////////////////////////////////////////////////////////////////
void LowerMultiplicityTree::beginJob() 
{
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"       ,&run_         ,"run_/I");
  outTree_->Branch("evtNo"       ,&evt_         ,"evt_/I");
  outTree_->Branch("nvtx"        ,&nVtx_        ,"nVtx_/I");
  outTree_->Branch("ht"          ,&ht_          ,"ht_/F");
  outTree_->Branch("pt"          ,&pt_          ,"pt_[8]/F");
  outTree_->Branch("eta"         ,&eta_         ,"eta_[8]/F");
  outTree_->Branch("phi"         ,&phi_         ,"phi_[8]/F");
  outTree_->Branch("njet"        ,&njet_        ,"njet_/I");
}
//////////////////////////////////////////////////////////////////////////////////////////
void LowerMultiplicityTree::endJob() 
{
  ;
}
//////////////////////////////////////////////////////////////////////////////////////////
void LowerMultiplicityTree::initialize()
{
  run_          = -999;
  evt_          = -999;
  nVtx_         = -999;
  ht_           = -999;
  njet_         = -999;

  for(int i=0;i<8;i++) {
    pt_[i]   = -999;
    eta_[i]  = -999;
    phi_[i]  = -999;
  }  
}
//////////////////////////////////////////////////////////////////////////////////////////

void LowerMultiplicityTree::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
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
  edm::Handle<GenJetCollection> genjets;
  edm::Handle<GenParticleCollection> genParticles;
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;

  initialize();

  int intpu(0);
  if (!iEvent.isRealData()) {
    if (srcPU_ != "") {
      iEvent.getByLabel(srcPU_,PupInfo); 
      for(vector<PileupSummaryInfo>::const_iterator ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu) {
        if (ipu->getBunchCrossing() == 0) {
          intpu += ipu->getPU_NumInteractions();
        } 
      }
    }
  } 
  bool cut_vtx = (recVtxs->size() > 0);
  bool cut_njets = (pat_jets.size() > 4);
  bool cut_pt8=false;
  
  if (cut_vtx && cut_njets) {
    bool cutID(true),cut_eta(true),cut_pt(true),cut_pu(true);
    vector<LorentzVector> P4(0); 
    LorentzVector P48J(0,0,0,0);
    int N(0);
    double HT(0.0);
    for(edm::View<pat::Jet>::const_iterator ijet = pat_jets.begin();ijet != pat_jets.end() && N < 9; ++ijet) { 
      bool id(false);
      double chf = ijet->chargedHadronEnergyFraction();
      double nhf = ijet->neutralHadronEnergyFraction() + ijet->HFHadronEnergyFraction();
      double phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double muf = ijet->muonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      int chm = ijet->chargedHadronMultiplicity();
      int npr = ijet->chargedMultiplicity() + ijet->neutralMultiplicity();
      id = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(ijet->eta())<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || fabs(ijet->eta())>2.4));
      double beta = ijet->userFloat("beta");
      if (N<8){
	cut_pu  *= (beta > betaMin_);
	cutID   *= id;
	cut_pt  *= (ijet->pt() > ptMin_);
	cut_eta *= (fabs(ijet->eta()) < etaMax_);
	HT += ijet->pt();
	pt_[N]  = ijet->pt();
	phi_[N] = ijet->phi();
	eta_[N] = ijet->eta();
	
	
      }
      
      if (ijet->pt() < ptMin_ && N==8){
	cut_pt=true;
      }
      N++;
      
      
    }// jet loop
   
    if (cutID && cut_pu && cut_eta && cut_pt) {
      ht_     = HT;
      run_    = iEvent.id().run();
      evt_    = iEvent.id().event();
      nVtx_   = recVtxs->size();
      if (pt_[4]>30. && pt_[5]<30.){
	njet_=5;
      }
      else if (pt_[5]>30. && pt_[6]<30.){
	njet_=6;
      }
      else if (pt_[6]>30. && pt_[7]<30.){
	njet_=7;
      }
      else if (pt_[7]>30. && cut_pt8){
	njet_=8;
      }    
      
      outTree_->Fill();
    }// if jet cuts
  }// if vtx and jet multi

}
//////////////////////////////////////////////////////////////////////////////////////////
LowerMultiplicityTree::~LowerMultiplicityTree() 
{
}

DEFINE_FWK_MODULE(LowerMultiplicityTree);
