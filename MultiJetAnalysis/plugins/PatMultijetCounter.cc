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

#include "KKousour/MultiJetAnalysis/plugins/PatMultijetCounter.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

using namespace std;
using namespace reco;

PatMultijetCounter::PatMultijetCounter(edm::ParameterSet const& cfg) 
{
  srcJets_ = cfg.getParameter<edm::InputTag>   ("jets");
  isMC_    = cfg.getUntrackedParameter<bool>   ("isMC",false);
  etaMax_  = cfg.getParameter<double>          ("etaMAX");
  ptMin_   = cfg.getParameter<double>          ("ptMIN");
  htMin_   = cfg.getParameter<double>          ("htMIN");
  betaMin_ = cfg.getParameter<double>          ("betaMIN");
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetCounter::beginJob() 
{
  hJetMulti_ = fs_->make<TH1F>("jetMulti","jetMulti",20,0,20);
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetCounter::endJob() 
{
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetCounter::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(srcJets_,jets);
  edm::View<pat::Jet> pat_jets = *jets;

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("goodOfflinePrimaryVertices",recVtxs);

  bool cut_vtx = (recVtxs->size() > 0);

  if (cut_vtx) {
    bool cut_eta(true),cut_pt(true),cut_pu(true);
    int N(0);
    double HT(0.0);
    for(edm::View<pat::Jet>::const_iterator ijet = pat_jets.begin();ijet != pat_jets.end(); ++ijet) { 
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
      cut_pu  = (beta > betaMin_);
      cut_pt  = (ijet->pt() > ptMin_);
      cut_eta = (fabs(ijet->eta()) < etaMax_);
      if (cut_pu && cut_pt && cut_eta && id) {
        HT += ijet->pt();
        N++;
      }
    }// jet loop
    if (HT > htMin_) {
      hJetMulti_->Fill(N);
    }
  }// if vtx
}
//////////////////////////////////////////////////////////////////////////////////////////
PatMultijetCounter::~PatMultijetCounter() 
{
}

DEFINE_FWK_MODULE(PatMultijetCounter);
