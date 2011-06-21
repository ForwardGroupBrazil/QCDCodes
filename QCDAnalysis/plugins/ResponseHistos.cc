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

#include "KKousour/QCDAnalysis/plugins/ResponseHistos.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "KKousour/QCDAnalysis/interface/QCDJet.h"
#include "KKousour/QCDAnalysis/interface/QCDEvent.h"
#include "KKousour/QCDAnalysis/interface/QCDEventHdr.h"
#include "KKousour/QCDAnalysis/interface/QCDCaloJet.h"
#include "KKousour/QCDAnalysis/interface/QCDPFJet.h"

using namespace std;

ResponseHistos::ResponseHistos(edm::ParameterSet const& cfg) 
{
  mPtBND     = cfg.getParameter<std::vector<double> > ("ptBnd");
  mEtaBND    = cfg.getParameter<std::vector<double> > ("etaBnd");
  mMaxDR     = cfg.getParameter<double> ("maxDR");
  mNEvents   = cfg.getParameter<int>  ("nEvents");
  mFileName  = cfg.getParameter<std::string> ("filename");
  mTreeName  = cfg.getParameter<std::string> ("treename");
  mDirName   = cfg.getParameter<std::string> ("dirname");
}
//////////////////////////////////////////////////////////////////////////////////////////
void ResponseHistos::beginJob() 
{
  char name[1000];
  //----------- Response Vs Eta Histos -----------
  for (unsigned ipt=0;ipt<mPtBND.size()-1;ipt++) {
    sprintf(name,"CaloRspVsEta_PtBin%d",ipt);
    mCaloRspVsEta[ipt] = fs->make<TH2F>(name,name,200,-5,5,500,0,2);
    mCaloRspVsEta[ipt]->Sumw2();
    sprintf(name,"PFRspVsEta_PtBin%d",ipt);
    mPFRspVsEta[ipt] = fs->make<TH2F>(name,name,200,-5,5,500,0,2);
    mPFRspVsEta[ipt]->Sumw2();
  }
  //----------- Response Vs Pt Histos ------------
  for (unsigned ieta=0;ieta<mEtaBND.size()-1;ieta++) {
    sprintf(name,"CaloRspVsPt_EtaBin%d",ieta);
    mCaloRspVsPt[ieta] = fs->make<TH2F>(name,name,140,0,3500,500,0,2);
    mCaloRspVsPt[ieta]->Sumw2();
    sprintf(name,"PFRspVsPt_EtaBin%d",ieta);
    mPFRspVsPt[ieta] = fs->make<TH2F>(name,name,140,0,3500,500,0,2);
    mPFRspVsPt[ieta]->Sumw2();
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void ResponseHistos::endJob() 
{
  mInf->Close();
}
//////////////////////////////////////////////////////////////////////////////////////////
void ResponseHistos::analyze(edm::Event const& event, edm::EventSetup const& iSetup) 
{
  mInf = TFile::Open(mFileName.c_str());
  mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());
  mTree = (TTree*)mDir->Get(mTreeName.c_str());
  mEvent = new QCDEvent();
  TBranch *branch = mTree->GetBranch("events");
  branch->SetAddress(&mEvent);
  
  unsigned NEntries = mTree->GetEntries();
  cout<<"File: "<<mFileName<<endl;
  cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
  int decade = 0;
  unsigned NN = NEntries;
  if (mNEvents > -1)
    NN = (unsigned)mNEvents;
  for(unsigned i=0;i<NN;i++) {
    double progress = 10.0*i/(1.0*NEntries);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;          
    mTree->GetEntry(i);
    double wt = mEvent->evtHdr().weight();
    //--------- Loop over PFJets --------------
    for(unsigned j=0;j<mEvent->nPFJets();j++) {
      double rGen   = (mEvent->pfjet(j)).genR();
      double ptGen  = (mEvent->pfjet(j)).genpt();
      double etaGen = (mEvent->pfjet(j)).geneta();
      if (rGen < mMaxDR) {
        int etaBin = getBin(etaGen,mEtaBND);
        int ptBin  = getBin(ptGen,mPtBND);
        if (etaBin > -1 && ptBin > -1) { 
          double rsp = (mEvent->pfjet(j)).ptCor()/ptGen;
          mPFRspVsEta[ptBin]->Fill(etaGen,rsp,wt);
          mPFRspVsPt[etaBin]->Fill(ptGen,rsp,wt);
        }
      }
    }
    //--------- Loop over CaloJets ------------
    for(unsigned j=0;j<mEvent->nCaloJets();j++) {
      double rGen   = (mEvent->calojet(j)).genR();
      double ptGen  = (mEvent->calojet(j)).genpt();
      double etaGen = (mEvent->calojet(j)).geneta();
      if (rGen < mMaxDR) {
        int etaBin = getBin(etaGen,mEtaBND);
        int ptBin  = getBin(ptGen,mPtBND);
        if (etaBin > -1 && ptBin > -1) {
          double rsp = (mEvent->calojet(j)).ptCor()/ptGen;
          mCaloRspVsEta[ptBin]->Fill(etaGen,rsp,wt);
          mCaloRspVsPt[etaBin]->Fill(ptGen,rsp,wt);
        }
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
int ResponseHistos::getBin(double x, const std::vector<double>& boundaries)
{
  int i;
  int n = boundaries.size()-1;
  if (x<boundaries[0] || x>=boundaries[n])
    return -1;
  for(i=0;i<n;i++)
   {
     if (x>=boundaries[i] && x<boundaries[i+1])
       return i;
   }
  return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////
ResponseHistos::~ResponseHistos() 
{
}

DEFINE_FWK_MODULE(ResponseHistos);
