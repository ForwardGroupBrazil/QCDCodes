#include <iostream>
#include "TString.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
using std::cin;
using std::cout;
using std::endl;
void FillHistograms(TString FileName, bool ApplyTriggerSel, bool ApplyPuWeight, bool isMC)
{
  cout<<"Opening files..............."<<endl;
  TFile *inf      = TFile::Open(FileName+".root");
  TFile *puf      = TFile::Open("/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/pileUpFile.root");
  TH1F *hPuWeight = (TH1F*)puf->Get("pileUpWeight");
  TFile *outf     = TFile::Open(FileName+"_histos.root","RECREATE");
  
  TTree *tr = (TTree*)inf->Get("Hbb/events");
  //---- define histograms ------------------
  cout<<"Booking histograms.........."<<endl;
  TH1F *hMbb          = new TH1F("hMbb","hMbb",100,0,1000);
  TH1F *hMbbCor       = new TH1F("hMbbCor","hMbbCor",100,0,1000);
  TH1F *hMqq          = new TH1F("hMqq","hMqq",100,0,4000);
  TH1F *hHtAll        = new TH1F("hHtAll","hHtAll",200,0,4000);
  TH1F *hdPhibb       = new TH1F("hdPhibb","hdPhibb",64,0,3.2001);
  TH1F *hEtaBoost     = new TH1F("hEtaBoost","hEtaBoost",60,-4,4);
  TH1F *hdEtaqq       = new TH1F("hdEtaqq","hdEtaqq",40,2,10);
  TH1F *hdEtaqqDiff   = new TH1F("hdEtaqqDiff","hdEtaqqDiff",60,0,3);
  TH1F *hsoftHt       = new TH1F("hsoftHt","hsoftHt",100,0,200);
  TH1F *hsoftMulti    = new TH1F("hsoftMulti","hsoftMulti",50,0,50);
  TH1F *hMet          = new TH1F("hMet","hMet",70,0,350);
  TH1F *hRho          = new TH1F("hRho","hRho",60,0,60);
  TH1F *hBDT          = new TH1F("hBDT","hBDT",100,-1,1);
  TH1F *hMLP          = new TH1F("hMLP","hMLP",68,-0.2,1.5);
  TH1F *hJetPt4Mqq    = new TH1F("hJetPt4Mqq","hJetPt4Mqq",50,0,0.25);
  TH1F *hMbbCut       = new TH1F("hMbbCut","hMbbCut",25,0,500);
  TH1F *hMbbCorCut    = new TH1F("hMbbCorCut","hMbbCorCut",25,0,500);
  TH1F *hMqqCut       = new TH1F("hMqqCut","hMqqCut",50,0,4000);
  TH1F *hHtAllCut     = new TH1F("hHtAllCut","hHtAllCut",50,0,4000);
  TH1F *hdPhibbCut    = new TH1F("hdPhibbCut","hdPhibbCut",25,0,3.14159);
  TH1F *hEtaBoostCut  = new TH1F("hEtaBoostCut","hEtaBoostCut",30,-4,4);
  TH1F *hdEtaqqCut    = new TH1F("hdEtaqqCut","hdEtaqqCut",35,3,10);
  TH1F *hdEtaqqDiffCut= new TH1F("hdEtaqqDiffCut","hdEtaqqDiffCut",60,0,3);
  TH1F *hsoftHtCut    = new TH1F("hsoftHtCut","hsoftHtCut",50,0,200);
  TH1F *hsoftMultiCut = new TH1F("hsoftMultiCut","hsoftMultiCut",50,0,50);
  TH1F *hMetCut       = new TH1F("hMetCut","hMetCut",35,0,350);
  TH1F *hBDTCut       = new TH1F("hBDTCut","hBDTCut",50,-1,1);
  TH1F *hMLPCut       = new TH1F("hMLPCut","hMLPCut",34,-0.2,1.5);
  TH1F *hJetPt4MqqCut = new TH1F("hJetPt4MqqCut","hJetPt4MqqCut",25,0,0.25);
  
  TH1F *hJetPt[5],*hJetPuMvaBtag[5],*hJetPtBtag[5],*hJetEta[5],*hJetEtaBtag[5],*hJetPhi[5],*hJetBtag[5],*hJetQGL[5],*hJetChf[5],
  *hJetNhf[5],*hJetPhf[5],*hJetMuf[5],*hJetElf[5],*hJetPtD[5];
  TH1F *hJetPtCut[5],*hJetPtBtagCut[5],*hJetEtaCut[5],*hJetEtaBtagCut[5],*hJetPhiCut[5],*hJetBtagCut[5],*hJetQGLCut[5],*hJetChfCut[5],
  *hJetNhfCut[5],*hJetPhfCut[5],*hJetMufCut[5],*hJetElfCut[5],*hJetPtDCut[5];
  char name[1000];
  for(int j=0;j<5;j++) {
    sprintf(name,"hJetPt%d",j);
    hJetPt[j] = new TH1F(name,name,200,0,1000);
    sprintf(name,"hJetPuMvaBtag%d",j);
    hJetPuMvaBtag[j] = new TH1F(name,name,100,-1,1);
    sprintf(name,"hJetPtBtag%d",j);
    hJetPtBtag[j] = new TH1F(name,name,200,0,1000);
    sprintf(name,"hJetEta%d",j);
    hJetEta[j] = new TH1F(name,name,50,-5,5);
    sprintf(name,"hJetEtaBtag%d",j);
    hJetEtaBtag[j] = new TH1F(name,name,50,-5,5);
    sprintf(name,"hJetPhi%d",j);
    hJetPhi[j] = new TH1F(name,name,50,-3.14159,3.14159);
    sprintf(name,"hJetBtag%d",j);
    hJetBtag[j] = new TH1F(name,name,50,0,1);
    sprintf(name,"hJetQGL%d",j);
    hJetQGL[j] = new TH1F(name,name,50,0,1);
    sprintf(name,"hJetPhf%d",j);
    hJetPhf[j] = new TH1F(name,name,50,0,1);
    sprintf(name,"hJetNhf%d",j);
    hJetNhf[j] = new TH1F(name,name,50,0,1);
    sprintf(name,"hJetElf%d",j);
    hJetElf[j] = new TH1F(name,name,50,0,1);
    sprintf(name,"hJetMuf%d",j);
    hJetMuf[j] = new TH1F(name,name,50,0,1);
    sprintf(name,"hJetChf%d",j);
    hJetChf[j] = new TH1F(name,name,50,0,1);
    sprintf(name,"hJetPtD%d",j);
    hJetPtD[j] = new TH1F(name,name,50,0,1);
    
    sprintf(name,"hJetPtCut%d",j);
    hJetPtCut[j] = new TH1F(name,name,200,0,1000);
    sprintf(name,"hJetPtBtagCut%d",j);
    hJetPtBtagCut[j] = new TH1F(name,name,200,0,1000);
    sprintf(name,"hJetEtaCut%d",j);
    hJetEtaCut[j] = new TH1F(name,name,25,-5,5);
    sprintf(name,"hJetEtaBtagCut%d",j);
    hJetEtaBtagCut[j] = new TH1F(name,name,25,-5,5);
    sprintf(name,"hJetPhiCut%d",j);
    hJetPhiCut[j] = new TH1F(name,name,25,-3.14159,3.14159);
    sprintf(name,"hJetBtagCut%d",j);
    hJetBtagCut[j] = new TH1F(name,name,50,0,1);
    sprintf(name,"hJetQGLCut%d",j);
    hJetQGLCut[j] = new TH1F(name,name,25,0,1);
    sprintf(name,"hJetPhfCut%d",j);
    hJetPhfCut[j] = new TH1F(name,name,50,0,1);
    sprintf(name,"hJetNhfCut%d",j);
    hJetNhfCut[j] = new TH1F(name,name,50,0,1);
    sprintf(name,"hJetElfCut%d",j);
    hJetElfCut[j] = new TH1F(name,name,50,0,1);
    sprintf(name,"hJetMufCut%d",j);
    hJetMufCut[j] = new TH1F(name,name,50,0,1);
    sprintf(name,"hJetChfCut%d",j);
    hJetChfCut[j] = new TH1F(name,name,50,0,1);
    sprintf(name,"hJetPtDCut%d",j);
    hJetPtDCut[j] = new TH1F(name,name,50,0,1);
  }
  cout<<"Defining variables.........."<<endl;
  float dEtaqq,dEtaqqEta,dPhibb,htAll,mqq,mbb,mbbCor,softHt,etaBoostqq;
  int btagIdx[5],softMulti,nvtx;
  vector<bool> *triggerResult(0);
  float jetQGL[5],jetBtag[5],jetPt[5],jetPuMva[5],jetPtD[5],jetEta[5],jetPhi[5],jetChf[5],jetNhf[5],jetPhf[5],jetElf[5],jetMuf[5];
  float met,rho;
  float BDT,MLP;
  
  tr->SetBranchAddress("jetQGLnew"     ,&jetQGL);
  tr->SetBranchAddress("btagIdx"       ,&btagIdx);
  tr->SetBranchAddress("mqq"           ,&mqq);
  tr->SetBranchAddress("mbb"           ,&mbb);
  tr->SetBranchAddress("mbbCor"        ,&mbbCor);
  tr->SetBranchAddress("dEtaqq"        ,&dEtaqq);
  tr->SetBranchAddress("dEtaqqEta"     ,&dEtaqqEta);
  tr->SetBranchAddress("dPhibb"        ,&dPhibb);
  tr->SetBranchAddress("etaBoostqq"    ,&etaBoostqq);
  tr->SetBranchAddress("htAll"         ,&htAll);
  tr->SetBranchAddress("nSoftTrackJets",&softMulti);
  tr->SetBranchAddress("softHt"        ,&softHt);
  tr->SetBranchAddress("rho"           ,&rho);
  tr->SetBranchAddress("nvtx"          ,&nvtx);
  tr->SetBranchAddress("met"           ,&met);
  tr->SetBranchAddress("BDT"           ,&BDT);
  tr->SetBranchAddress("MLP"           ,&MLP);
  tr->SetBranchAddress("jetPt"         ,&jetPt);
  tr->SetBranchAddress("jetPuMva"      ,&jetPuMva); 
  tr->SetBranchAddress("jetPtD"        ,&jetPtD);
  tr->SetBranchAddress("jetEta"        ,&jetEta);
  tr->SetBranchAddress("jetPhi"        ,&jetPhi);
  tr->SetBranchAddress("jetBtag"       ,&jetBtag);
  tr->SetBranchAddress("jetChf"        ,&jetChf);
  tr->SetBranchAddress("jetNhf"        ,&jetNhf);
  tr->SetBranchAddress("jetPhf"        ,&jetPhf);
  tr->SetBranchAddress("jetElf"        ,&jetElf);
  tr->SetBranchAddress("jetMuf"        ,&jetMuf);
  tr->SetBranchAddress("triggerResult" ,&triggerResult);
  
  int decade(0);
  int NN = tr->GetEntries();
  cout<<"Reading "<<NN<<" entries"<<endl;
  for(int i=0;i<tr->GetEntries();i++) {
    double progress = 10.0*i/(1.0*NN);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;
    tr->GetEntry(i);
    bool cut_trigger(true);
    bool cut_btag = ((jetBtag[btagIdx[0]] > 0) && (jetBtag[btagIdx[1]] > 0));
    if (ApplyTriggerSel) {
      if (isMC) {
        cut_trigger = ((*triggerResult)[5] || (*triggerResult)[7]);
      }
      else {
        cut_trigger = ((*triggerResult)[0] || (*triggerResult)[1]);
      }
    }
    if (!cut_trigger) continue;
    if (!cut_btag) continue;
    float wt(1.0);
    float wtPu(1.0);
    if (ApplyPuWeight) {
      int bin = hPuWeight->FindBin(rho);
      wtPu = hPuWeight->GetBinContent(bin);
    }
    wt *= wtPu;
    if (dPhibb > 2.0) continue;
    if (softHt > 400) continue;
    bool PuId(true);
    for(int j=0;j<4;j++) {
      float eta = fabs(jetEta[j]);
      if (eta < 2.5) {
        PuId *= (jetPuMva[j] > -0.38);
      }
      else if (eta >= 2.5 && eta < 2.75) {
        PuId *= (jetPuMva[j] > -0.32);
      }
      else if (eta >= 2.75 && eta < 3.0) {
        PuId *= (jetPuMva[j] > -0.14);
      }
      else {
        PuId *= (jetPuMva[j] > -0.48);
      }
    }
    if (!PuId) continue;
    hMbb->Fill(mbb,wt);
    hMbbCor->Fill(mbbCor,wt);
    hMqq->Fill(mqq,wt);
    hdEtaqq->Fill(dEtaqq,wt);
    hdEtaqqDiff->Fill(dEtaqqEta-dEtaqq,wt);
    hdPhibb->Fill(dPhibb,wt);
    hHtAll->Fill(htAll,wt);
    hEtaBoost->Fill(etaBoostqq,wt);
    hsoftHt->Fill(softHt,wt);
    hsoftMulti->Fill(softMulti,wt);
    hMet->Fill(met,wt);
    hMLP->Fill(MLP,wt);
    hBDT->Fill(BDT,wt);
    hRho->Fill(rho,wt);
    if (jetPt[4] > 0) {
      hJetPt4Mqq->Fill(jetPt[4]/mqq,wt);
    }
    for(int j=0;j<5;j++) {
      if (jetPt[j]<0) continue;
      int ib = btagIdx[j];
      hJetPt[j]->Fill(jetPt[j],wt);
      hJetPuMvaBtag[j]->Fill(jetPuMva[ib],wt);
      hJetPtBtag[j]->Fill(jetPt[ib],wt);
      hJetEta[j]->Fill(jetEta[j],wt);
      hJetEtaBtag[j]->Fill(jetEta[ib],wt);
      hJetPhi[j]->Fill(jetPhi[j],wt);
      hJetChf[j]->Fill(jetChf[j],wt);
      hJetNhf[j]->Fill(jetNhf[j],wt);
      hJetPhf[j]->Fill(jetPhf[j],wt);
      hJetMuf[j]->Fill(jetMuf[j],wt);
      hJetElf[j]->Fill(jetElf[j],wt);
      hJetPtD[j]->Fill(jetPtD[j],wt);
      hJetBtag[j]->Fill(jetBtag[ib],wt);
      hJetQGL[j]->Fill(jetQGL[ib],wt);
    }
    if (MLP > 0.9 && dPhibb < 2.0) {
      hMbbCut->Fill(mbb,wt);
      hMbbCorCut->Fill(mbbCor,wt);
      hMqqCut->Fill(mqq,wt);
      hdEtaqqCut->Fill(dEtaqq,wt);
      hdEtaqqDiffCut->Fill(dEtaqqEta-dEtaqq,wt);
      hdPhibbCut->Fill(dPhibb,wt);
      hHtAllCut->Fill(htAll,wt);
      hEtaBoostCut->Fill(etaBoostqq,wt);
      hsoftHtCut->Fill(softHt,wt);
      hsoftMultiCut->Fill(softMulti,wt);
      hMetCut->Fill(met,wt);
      hMLPCut->Fill(MLP,wt);
      hBDTCut->Fill(BDT,wt);
      if (jetPt[4] > 0) {
        hJetPt4MqqCut->Fill(jetPt[4]/mqq,wt);
      }
      for(int j=0;j<5;j++) {
        if (jetPt[j]<0) continue;
        int ib = btagIdx[j];
        hJetPtCut[j]->Fill(jetPt[j],wt);
        hJetPtBtagCut[j]->Fill(jetPt[ib],wt);
        hJetEtaCut[j]->Fill(jetEta[j],wt);
        hJetEtaBtagCut[j]->Fill(jetEta[ib],wt);
        hJetPhiCut[j]->Fill(jetPhi[j],wt);
        hJetChfCut[j]->Fill(jetChf[j],wt);
        hJetNhfCut[j]->Fill(jetNhf[j],wt);
        hJetPhfCut[j]->Fill(jetPhf[j],wt);
        hJetMufCut[j]->Fill(jetMuf[j],wt);
        hJetElfCut[j]->Fill(jetElf[j],wt);
        hJetPtDCut[j]->Fill(jetPtD[j],wt);
        hJetBtagCut[j]->Fill(jetBtag[ib],wt);
        hJetQGLCut[j]->Fill(jetQGL[ib],wt);
      }
    }// cut
  }// tree loop  
  outf->Write();
  inf->Close();
  outf->Close();
}





