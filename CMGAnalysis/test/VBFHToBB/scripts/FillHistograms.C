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
void FillHistograms(TString FileName, bool ApplyTriggerSel, bool isMC)
{
  cout<<"Opening files..............."<<endl;
  TFile *inf      = TFile::Open(FileName+".root");
  //TFile *puf      = TFile::Open("/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/pileUpFile.root");
  //TH1F *hPuWeight = (TH1F*)puf->Get("pileUpWeight");
  TFile *outf     = TFile::Open(FileName+"_histos.root","RECREATE");
  
  TTree *tr = (TTree*)inf->Get("Hbb/events");
  //---- define histograms ------------------
  const int NVAR = 15;
  TString var[NVAR] = {"mbb","mbbCor","mqq","dPhibb","etaBoostqq","dEtaqq","dEtaqqDiff","softHt","softMulti",
                       "cosTheta","cosAlpha","met","rho","MLP","nVtx"}; 
  int NBINS[NVAR]   = {100,100,100,80,60,40,60,100,50,50,50,70,60,100,50};
  double XMIN[NVAR] = {0,0,0,0,-4,2,0,0,0,-1,-1,0,0,-0.5,0};
  double XMAX[NVAR] = {1000,1000,4000,3.2,4,10,3,200,50,1.0001,1.0001,350,60,1.5,50};
  TH1F *hVar[NVAR],*hVarCut[NVAR];
  cout<<"Booking histograms.........."<<endl;
  for(int ivar=0;ivar<NVAR;ivar++) {
    hVar[ivar]    = new TH1F("h"+var[ivar],"h"+var[ivar],NBINS[ivar],XMIN[ivar],XMAX[ivar]);
    hVarCut[ivar] = new TH1F("h"+var[ivar]+"Cut","h"+var[ivar]+"Cut",NBINS[ivar]/4,XMIN[ivar],XMAX[ivar]);
  }
  cout<<"Booking jet histograms.........."<<endl;
  const int NJETVAR = 14;
  TString varJet[NJETVAR] = {"jetPt","jetPtBtag","jetEta","jetEtaBtag","jetPhi","jetBtag","jetQGL","jetChf",
                       "jetNhf","jetPhf","jetMuf","jetElf","jetPtD","jetPuMva"}; 
  int NJETBINS[NJETVAR]   = {200,200,50,50,50,50,50,50,50,50,50,50,50,100};
  double XMINJET[NJETVAR] = {0,0,-5,-5,-3.14159,0,0,0,0,0,0,0,0,-1};
  double XMAXJET[NJETVAR] = {1000,1000,5,5,3.14159,1,1,1,1,1,1,1,1,1};
  TH1F *hJetVar[NJETVAR][5],*hJetVarCut[NJETVAR][5];
  char name[1000];
  for(int ivar=0;ivar<NJETVAR;ivar++) {
    for(int j=0;j<5;j++) {
      sprintf(name,"h%s%d",varJet[ivar].Data(),j);
      hJetVar[ivar][j]    = new TH1F(name,name,NJETBINS[ivar],XMINJET[ivar],XMAXJET[ivar]);
      sprintf(name,"h%sCut%d",varJet[ivar].Data(),j);
      hJetVarCut[ivar][j] = new TH1F(name,name,NJETBINS[ivar]/4,XMINJET[ivar],XMAXJET[ivar]);
    }
  }
  cout<<"Defining variables.........."<<endl;
  float dEtaqq,dEtaqqEta,dPhibb,htAll,mqq,mbb,mbbCor,softHt,etaBoostqq,puWt,cosTheta,cosAlpha;
  int btagIdx[5],softMulti,nvtx,puId[3];
  vector<bool> *triggerResult(0);
  float jetQGL[5],jetBtag[5],jetPt[5],jetPuMva[5],jetPtD[5],jetEta[5],jetPhi[5],jetChf[5],jetNhf[5],jetPhf[5],jetElf[5],jetMuf[5];
  float met,rho;
  float MLP;
  
  tr->SetBranchAddress("jetQGLnew"     ,&jetQGL);
  tr->SetBranchAddress("btagIdx"       ,&btagIdx);
  tr->SetBranchAddress("puId"          ,&puId);
  tr->SetBranchAddress("puWt"          ,&puWt);
  tr->SetBranchAddress("mqq"           ,&mqq);
  tr->SetBranchAddress("mbb"           ,&mbb);
  tr->SetBranchAddress("mbbCor"        ,&mbbCor);
  tr->SetBranchAddress("cosTheta"      ,&cosTheta);
  tr->SetBranchAddress("cosAlpha"      ,&cosAlpha);
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
    if (!cut_btag)    continue;
    if (dPhibb > 2.0) continue;
    if (softHt > 400) continue;
    if (puId[0] == 0) continue;
    float wt = puWt;
    double x[NVAR] = {mbb,mbbCor,mqq,dPhibb,etaBoostqq,dEtaqq,dEtaqqEta-dEtaqq,softHt,softMulti,cosTheta,cosAlpha,
                 met,rho,MLP,nvtx};
    for(int ivar=0;ivar<NVAR;ivar++) {
      hVar[ivar]->Fill(x[ivar],wt);
      if (MLP > 0.8) {
        hVarCut[ivar]->Fill(x[ivar],wt);
      }
    }
    for(int j=0;j<5;j++) {
      if (jetPt[j]<0) continue;
      int ib = btagIdx[j]; 
      double xJet[NJETVAR] = {jetPt[j],jetPt[ib],jetEta[j],jetEta[ib],jetPhi[j],jetBtag[ib],jetQGL[ib],jetChf[j],jetNhf[j],jetPhf[j],
                              jetMuf[j],jetElf[j],jetPtD[j],jetPuMva[ib]}; 
      for(int ivar=0;ivar<NJETVAR;ivar++) {
        hJetVar[ivar][j]->Fill(xJet[ivar],wt);
        if (MLP > 0.8) {
          hJetVarCut[ivar][j]->Fill(xJet[ivar],wt);
        }
      }
    }
  }// tree loop  
  outf->Write();
  inf->Close();
  outf->Close();
}





