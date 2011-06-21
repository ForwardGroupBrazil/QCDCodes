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

#include "KKousour/QCDAnalysis/plugins/DijetGENSearchHistos.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"


using namespace std;

DijetGENSearchHistos::DijetGENSearchHistos(edm::ParameterSet const& cfg) 
{
  mMassBND   = cfg.getParameter<std::vector<double> > ("massBnd");
  mMinMass   = cfg.getParameter<double> ("minMass");
  mMinPt     = cfg.getParameter<double> ("minPt");
  mMaxEta    = cfg.getParameter<double> ("maxEta");
  mMaxDeta   = cfg.getParameter<double> ("maxDeta");
  mChiIN     = cfg.getParameter<double> ("chiIN");
  mChiOUT     = cfg.getParameter<double> ("chiOUT");
  mGenJetsName = cfg.getParameter<std::string> ("genjets");
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetGENSearchHistos::beginJob() 
{
  char name[1000];
  double auxMass[500];
  for(unsigned im=0;im<mMassBND.size();im++) {
    auxMass[im] = mMassBND[im];
  }
  sprintf(name,"JetPt");
  mhPt = fs->make<TH1F>(name,name,3500,0,3500);
  mhPt->Sumw2();
  sprintf(name,"JetPtJJ");
  mhPtJJ = fs->make<TH1F>(name,name,3500,0,3500);
  mhPtJJ->Sumw2();
  sprintf(name,"JetEta");
  mhEta = fs->make<TH1F>(name,name,300,-3,3);
  mhEta->Sumw2();
  sprintf(name,"Deta");
  mhDeta = fs->make<TH1F>(name,name,300,-3,3);
  mhDeta->Sumw2();
  sprintf(name,"etaBoost");
  mhEtaBoost = fs->make<TH1F>(name,name,300,-3,3);
  mhEtaBoost->Sumw2();
  sprintf(name,"Dphi");
  mhDphi = fs->make<TH1F>(name,name,200,0,3.15);
  mhDphi->Sumw2();
  sprintf(name,"Mass");
  mhMass = fs->make<TH1F>(name,name,mMassBND.size()-1,auxMass);
  mhMass->Sumw2();
  sprintf(name,"MassIN");
  mhMassIN = fs->make<TH1F>(name,name,mMassBND.size()-1,auxMass);
  mhMassIN->Sumw2();
  sprintf(name,"MassOUT");
  mhMassOUT = fs->make<TH1F>(name,name,mMassBND.size()-1,auxMass);
  mhMassOUT->Sumw2();
  sprintf(name,"Chi");
  mhChi = fs->make<TH1F>(name,name,100,0,20);
  mhChi->Sumw2();
  sprintf(name,"ChiIN");
  mhChiIN = fs->make<TH1F>(name,name,100,0,20);
  mhChiIN->Sumw2();
  sprintf(name,"ChiOUT");
  mhChiOUT = fs->make<TH1F>(name,name,100,0,20);
  mhChiOUT->Sumw2();
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetGENSearchHistos::endJob() 
{
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetGENSearchHistos::analyze(edm::Event const& event, edm::EventSetup const& iSetup) 
{ 
  edm::Handle<reco::GenJetCollection> genjets;
  event.getByLabel(mGenJetsName,genjets);
    
  if (genjets->size() > 1) {
    double deta = (*genjets)[0].eta()-(*genjets)[1].eta();
    double etaBoost = 0.5*((*genjets)[0].eta()+(*genjets)[1].eta());
    double dphi = (*genjets)[0].phi()-(*genjets)[1].phi();
    double mass = ((*genjets)[0].p4()+(*genjets)[1].p4()).mass();
    double ptJJ = ((*genjets)[0].p4()+(*genjets)[1].p4()).pt();
    double cosThetaStar = tanh(0.5*deta);
    double chi = (1+fabs(cosThetaStar))/(1-fabs(cosThetaStar));
    bool cutPt   = ((*genjets)[0].pt() >= mMinPt && (*genjets)[1].pt() >= mMinPt);
    bool cutEta  = ((*genjets)[0].eta() <= mMaxEta && (*genjets)[1].eta() <= mMaxEta);
    bool cutDeta = (fabs(deta) <= mMaxDeta);
    bool cutMass = (mass >= mMinMass);
    bool cutIN = (chi <= mChiIN);
    bool cutOUT = (chi > mChiIN && chi <= mChiOUT);
    if (cutPt && cutEta && cutMass) {
      mhPtJJ->Fill(ptJJ);
      mhChi->Fill(chi);
      mhEtaBoost->Fill(etaBoost); 
      if (cutIN) {
        mhMassIN->Fill(mass);
        mhChiIN->Fill(chi);
      }
      if (cutOUT) {
        mhMassOUT->Fill(mass);
        mhChiOUT->Fill(chi);
      }
      if (cutDeta) {
        mhMass->Fill(mass);
        mhDeta->Fill(deta);
        mhDphi->Fill(fabs(dphi));
        for(unsigned j=0;j<2;j++) {
          mhPt->Fill((*genjets)[j].pt());
          mhEta->Fill((*genjets)[j].eta()); 
        }
      }
    }
  } 
}
//////////////////////////////////////////////////////////////////////////////////////////
DijetGENSearchHistos::~DijetGENSearchHistos() 
{
}

DEFINE_FWK_MODULE(DijetGENSearchHistos);
