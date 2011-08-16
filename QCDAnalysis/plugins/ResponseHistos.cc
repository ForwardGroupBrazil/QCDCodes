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
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"


using namespace std;

ResponseHistos::ResponseHistos(edm::ParameterSet const& cfg) 
{
  mPtBND     = cfg.getParameter<std::vector<double> > ("ptBnd");
  mYBND      = cfg.getParameter<std::vector<double> > ("yBnd");
  mFineYBND  = cfg.getParameter<std::vector<double> > ("yFineBnd");
  mMaxDR     = cfg.getParameter<double>               ("maxDR");
  mNEvents   = cfg.getParameter<int>                  ("nEvents");
  mMaxJets   = cfg.getParameter<int>                  ("nJets");
  mFileName  = cfg.getParameter<std::string>          ("filename");
  mTreeName  = cfg.getParameter<std::string>          ("treename");
  mDirName   = cfg.getParameter<std::string>          ("dirname");
  mPFBiasCorName = cfg.getParameter<std::string>      ("PFBiasCorName");
}
//////////////////////////////////////////////////////////////////////////////////////////
void ResponseHistos::beginJob() 
{
  char name[1000];
  mGenYVsPt = fs->make<TH2F>("GenYVsPt","GenYVsPt",140,0,3500,200,-5,5);
  double vy[100];
  for(unsigned iy=0;iy<mFineYBND.size();iy++)
    vy[iy] = mFineYBND[iy];
  //----------- Response Vs Eta Histos -----------
  for (unsigned ipt=0;ipt<mPtBND.size()-1;ipt++) {
    sprintf(name,"CaloRspVsY_PtBin%d",ipt);
    mCaloRspVsY[ipt] = fs->make<TH2F>(name,name,200,-5,5,500,0,2);
    mCaloRspVsY[ipt]->Sumw2();
    sprintf(name,"PFRspVsY_PtBin%d",ipt);
    mPFRspVsY[ipt] = fs->make<TH2F>(name,name,200,-5,5,500,0,2);
    mPFRspVsY[ipt]->Sumw2();
    sprintf(name,"CaloDYVsY_PtBin%d",ipt);
    mCaloDYVsY[ipt] = fs->make<TProfile>(name,name,200,-5,5,-0.5,0.5);
    mCaloDYVsY[ipt]->Sumw2();
    sprintf(name,"PFDYVsY_PtBin%d",ipt);
    mPFDYVsY[ipt] = fs->make<TProfile>(name,name,200,-5,5,-0.5,0.5);
    mPFDYVsY[ipt]->Sumw2();
    sprintf(name,"CaloDYVsAbsY_PtBin%d",ipt);
    mCaloDYVsAbsY[ipt] = fs->make<TProfile>(name,name,mFineYBND.size()-1,vy,-0.5,0.5);
    mCaloDYVsAbsY[ipt]->Sumw2();
    sprintf(name,"PFDYVsAbsY_PtBin%d",ipt);
    mPFDYVsAbsY[ipt] = fs->make<TProfile>(name,name,mFineYBND.size()-1,vy,-0.5,0.5);
    mPFDYVsAbsY[ipt]->Sumw2();
    sprintf(name,"CaloPtVsAbsY_PtBin%d",ipt);
    mCaloPtVsAbsY[ipt] = fs->make<TProfile>(name,name,mFineYBND.size()-1,vy,0,3500);
    mCaloPtVsAbsY[ipt]->Sumw2();
    sprintf(name,"PFPtVsAbsY_PtBin%d",ipt);
    mPFPtVsAbsY[ipt] = fs->make<TProfile>(name,name,mFineYBND.size()-1,vy,0,3500);
    mPFPtVsAbsY[ipt]->Sumw2();
    if (mPFBiasCorName != "") {
      sprintf(name,"PFDCorYVsY_PtBin%d",ipt);
      mPFDCorYVsY[ipt] = fs->make<TProfile>(name,name,200,-5,5,-0.5,0.5);
      mPFDCorYVsY[ipt]->Sumw2();
    }
  }
  //----------- Response Vs Pt Histos ------------
  for (unsigned iy=0;iy<mYBND.size()-1;iy++) {
    sprintf(name,"CaloRspVsPt_YBin%d",iy);
    mCaloRspVsPt[iy] = fs->make<TH2F>(name,name,140,0,3500,500,0,2);
    mCaloRspVsPt[iy]->Sumw2();
    sprintf(name,"PFRspVsPt_YBin%d",iy);
    mPFRspVsPt[iy] = fs->make<TH2F>(name,name,140,0,3500,500,0,2);
    mPFRspVsPt[iy]->Sumw2();
    sprintf(name,"CaloDYVsPt_Ybin%d",iy);
    mCaloDYVsPt[iy] = fs->make<TH2F>(name,name,140,0,3500,600,-0.3,0.3);
    mCaloDYVsPt[iy]->Sumw2();
    sprintf(name,"PFDYVsPt_Ybin%d",iy);
    mPFDYVsPt[iy] = fs->make<TH2F>(name,name,140,0,3500,600,-0.3,0.3);
    mPFDYVsPt[iy]->Sumw2();
  }
  cout<<"Histograms booked"<<endl;
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
  JetCorrectorParameters *yBiasPar(0);
  JetCorrectionUncertainty *yBiasCor(0);
  if (mPFBiasCorName != "") {
    yBiasPar = new JetCorrectorParameters(mPFBiasCorName);
    yBiasCor = new JetCorrectionUncertainty(*yBiasPar);
  }
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
    //--------- Loop over GenJets --------------
    for(unsigned j=0;j<mEvent->nGenJets();j++) {
      mGenYVsPt->Fill(mEvent->genjet(j).pt(),mEvent->genjet(j).Rapidity(),wt);
    }
    //--------- Loop over PFJets --------------
    for(unsigned j=0;j<mEvent->nPFJets();j++) {
      double rGen   = (mEvent->pfjet(j)).genR();
      double ptGen  = (mEvent->pfjet(j)).genpt();
      double yGen   = (mEvent->pfjet(j)).genp4().Rapidity();
      if (rGen < mMaxDR) {
        int yBin = getBin(fabs(yGen),mYBND);
        int ptBin  = getBin(ptGen,mPtBND);
        if (ptBin > -1) {
          double rsp = (mEvent->pfjet(j)).ptCor()/ptGen;
          mPFDYVsY[ptBin]->Fill(yGen,mEvent->pfjet(j).y()-yGen,wt);
          mPFPtVsAbsY[ptBin]->Fill(fabs(yGen),mEvent->pfjet(j).ptCor(),wt); 
          mPFRspVsY[ptBin]->Fill(yGen,rsp,wt);
          if (mPFBiasCorName != "") {
            yBiasCor->setJetEta(mEvent->pfjet(j).y());
            yBiasCor->setJetPt(mEvent->pfjet(j).ptCor());
            mPFDCorYVsY[ptBin]->Fill(yGen,mEvent->pfjet(j).y()+yBiasCor->getUncertainty(true)-yGen,wt);
          }
          if (yGen < 0) {
            mPFDYVsAbsY[ptBin]->Fill(fabs(yGen),-(mEvent->pfjet(j).y()-yGen),wt);
          }
          else {
            mPFDYVsAbsY[ptBin]->Fill(fabs(yGen),mEvent->pfjet(j).y()-yGen,wt);  
          }
          if (yBin > -1) { 
            mPFDYVsPt[yBin]->Fill(ptGen,mEvent->pfjet(j).y()-yGen,wt);
            mPFRspVsPt[yBin]->Fill(ptGen,rsp,wt);
          }
        }
      }
    }
    //--------- Loop over CaloJets ------------
    for(unsigned j=0;j<mEvent->nCaloJets();j++) {
      double rGen   = (mEvent->calojet(j)).genR();
      double ptGen  = (mEvent->calojet(j)).genpt();
      double yGen   = (mEvent->calojet(j)).genp4().Rapidity();
      if (rGen < mMaxDR) {
        int yBin = getBin(fabs(yGen),mYBND);
        int ptBin  = getBin(ptGen,mPtBND);
        if (ptBin > -1) {
          double rsp = (mEvent->calojet(j)).ptCor()/ptGen;
          mCaloDYVsY[ptBin]->Fill(yGen,mEvent->calojet(j).y()-yGen,wt);
          mCaloPtVsAbsY[ptBin]->Fill(fabs(yGen),mEvent->calojet(j).ptCor(),wt);
          mCaloRspVsY[ptBin]->Fill(yGen,rsp,wt);
          if (yGen < 0)
            mCaloDYVsAbsY[ptBin]->Fill(fabs(yGen),-(mEvent->calojet(j).y()-yGen),wt);
          else
            mCaloDYVsAbsY[ptBin]->Fill(fabs(yGen),mEvent->calojet(j).y()-yGen,wt);
          if (yBin > -1) { 
            mCaloDYVsPt[yBin]->Fill(ptGen,mEvent->calojet(j).y()-yGen,wt);
            mCaloRspVsPt[yBin]->Fill(ptGen,rsp,wt);
          }
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
