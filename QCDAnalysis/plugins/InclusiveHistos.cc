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

#include "KKousour/QCDAnalysis/plugins/InclusiveHistos.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

InclusiveHistos::InclusiveHistos(edm::ParameterSet const& cfg) 
{
  mYBND            = cfg.getParameter<std::vector<double> >      ("yBnd");
  mPTBND           = cfg.getParameter<std::vector<double> >      ("ptBnd");
  mMinPt           = cfg.getParameter<std::vector<double> >      ("minPt");
  mMaxBetaStar     = cfg.getParameter<double>                    ("maxBetaStar");
  mMaxMETovSumET   = cfg.getParameter<double>                    ("maxMETovSumET"); 
  mFileNames       = cfg.getParameter<std::vector<std::string> > ("filenames");
  mTreeName        = cfg.getParameter<std::string>               ("treename");
  mDirName         = cfg.getParameter<std::string>               ("dirname");
  mLogName         = cfg.getParameter<std::string>               ("logname");
  mPUFileName      = cfg.getUntrackedParameter<std::string>      ("puFile","");
  mPUHistName      = cfg.getUntrackedParameter<std::string>      ("puHisto","");
  mTriggers        = cfg.getParameter<std::vector<std::string> > ("triggers");
  mIsMC            = cfg.getParameter<bool>                      ("isMC");
  mApplyHBEHfilter = cfg.getParameter<bool>                      ("applyHBEHfilter");
  mNEvents         = cfg.getParameter<int>                       ("nEvents"); 
  if (mTriggers.size() > 0) {
    if (mMinPt.size() != mTriggers.size()) {
      throw cms::Exception("InclusiveHistos: ")<<" Number of pt thresholds must be equal to the number of triggers\n";
    }
  }
  else {
    if (mMinPt.size() != 1) {
      throw cms::Exception("InclusiveHistos: ")<<" Number of pt thresholds must be equal to 1 when no triggers are defined\n";
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void InclusiveHistos::beginJob() 
{
  if (mIsMC && mPUFileName != "" && mPUHistName != "") {
    mPUf = TFile::Open(mPUFileName.c_str());
    mPUh = (TH1F*)mPUf->Get(mPUHistName.c_str());
  }
  TFile *mInf = TFile::Open(mFileNames[0].c_str());
  TDirectoryFile *mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());
  mEvent = new QCDEvent();
  char name[1000];
  //--------- trigger mapping -------------------
  TH1F *hTrigNames = (TH1F*)mDir->Get("TriggerNames");
  cout<<"Number of triggers = "<<mTriggers.size()<<endl;
  if (mTriggers.size() == 0)
    cout<<"No triggers set"<<endl;
  else {
    cout<<"Finding trigger mapping: "<<endl;
    mTrigIndex.clear();
    vector<int> versionTrig;
    for(unsigned itrig=0;itrig<mTriggers.size();itrig++) {
      versionTrig.clear();
      for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
        string ss = hTrigNames->GetXaxis()->GetBinLabel(ibin+1);
        //cout<<mTriggers[itrig]<<" "<<ss<<" "<<ss.find(mTriggers[itrig])<<endl;
        if (ss.find(mTriggers[itrig]) == 0) {
          versionTrig.push_back(ibin);
        }
      }
      if (versionTrig.size() > 0)
        mTrigIndex.push_back(versionTrig); 
      else {
        throw cms::Exception("InclusiveHistos: ")<<"The requested trigger ("<<mTriggers[itrig]<<") is not found ";
      }
    }
    for(unsigned itrig=0;itrig<mTriggers.size();itrig++) {
      cout<<mTriggers[itrig]<<" --> ";
      for(unsigned iv=0;iv<mTrigIndex[itrig].size();iv++) {
        cout<<mTrigIndex[itrig][iv]<<", ";
      }
      cout<<endl;
    }
  } 
  //--------- book histos ---------------------------------
  TFileDirectory mPFDir   = fs->mkdir(mDirName+"pf");
  TFileDirectory mCaloDir = fs->mkdir(mDirName+"calo");
  double auxPt[500],auxX[500];
  for(unsigned ipt=0;ipt<mPTBND.size();ipt++) {
    auxPt[ipt] = mPTBND[ipt];
    auxX[ipt] = mPTBND[ipt]/1000.; 
  }
  if (mIsMC) {
    mGenPU = fs->make<TH1F>("GenPU","GenPU",50,0,50);
  }
  mNPV = fs->make<TH1F>("NPV","NPV",50,0,50);
  mPVx = fs->make<TH1F>("PVx","PVx",200,-0.1,0.1);
  mPVy = fs->make<TH1F>("PVy","PVy",200,-0.1,0.1);
  mPVz = fs->make<TH1F>("PVz","PVz",200,-30,30);
  mBSx = fs->make<TH1F>("BSx","BSx",200,-0.1,0.1);
  mBSy = fs->make<TH1F>("BSy","BSy",200,-0.1,0.1);
  mBSz = fs->make<TH1F>("BSz","BSz",200,-0.2,0.2);
  mNvtxVsRun = fs->make<TProfile>("nVtxVsRun","nVtxVsRun",3,0,3,0,100);
  mNvtxVsRun->SetBit(TH1::kCanRebin);
  mCaloRhoVsRun = fs->make<TProfile>("CaloRhoVsRun","CaloRhoVsRun",3,0,3,0,50);
  mCaloRhoVsRun->SetBit(TH1::kCanRebin);
  mPFRhoVsRun = fs->make<TProfile>("PFRhoVsRun","PFRhoVsRun",3,0,3,0,50);
  mPFRhoVsRun->SetBit(TH1::kCanRebin);
  mCaloRhoVsNPV = fs->make<TProfile>("CaloRhoVsNPV","CaloRhoVsNPV",50,0,50,0,50);
  mPFRhoVsNPV = fs->make<TProfile>("PFRhoVsNPV","PFRhoVsNPV",50,0,50,0,50);
  int Ntrig = (int)mTriggers.size();
  for(int itrig=0;itrig<TMath::Max(Ntrig,1);itrig++) {
    string ss("");
    if (Ntrig > 0) {
      ss = "_"+mTriggers[itrig].erase(mTriggers[itrig].find("v")-1,mTriggers[itrig].find("v"));
    }
    sprintf(name,"METovSUMET_%s",ss.c_str());
    mPFMETovSUMET[itrig]   = mPFDir.make<TH1F>(name,name,100,0,1.0001);
    mCaloMETovSUMET[itrig] = mCaloDir.make<TH1F>(name,name,100,0,1.0001);
    sprintf(name,"JetMultiplicity%s",ss.c_str());
    mPFJetMulti[itrig]   = mPFDir.make<TH1F>(name,name,20,0,20);
    mCaloJetMulti[itrig] = mCaloDir.make<TH1F>(name,name,20,0,20);
    for(unsigned iy=0;iy<mYBND.size()-1;iy++) {
      //-------- Njets histograms -------------------------
      sprintf(name,"NPFJetsVsRun_Ybin%d%s",iy,ss.c_str());
      mNPFJets[itrig][iy] = mPFDir.make<TH1F>(name,name,3,0,3);
      mNPFJets[itrig][iy]->SetBit(TH1::kCanRebin);
      mNPFJets[itrig][iy]->Sumw2();
      sprintf(name,"NPFNormJetsVsRun_Ybin%d%s",iy,ss.c_str());
      mNPFNormJets[itrig][iy] = mPFDir.make<TH1F>(name,name,3,0,3);
      mNPFNormJets[itrig][iy]->SetBit(TH1::kCanRebin);
      mNPFNormJets[itrig][iy]->Sumw2();
      sprintf(name,"NCaloJetsVsRun_Ybin%d%s",iy,ss.c_str());
      mNCaloJets[itrig][iy] = mCaloDir.make<TH1F>(name,name,3,0,3);
      mNCaloJets[itrig][iy]->SetBit(TH1::kCanRebin);
      mNCaloJets[itrig][iy]->Sumw2();
      sprintf(name,"NCaloNormJetsVsRun_Ybin%d%s",iy,ss.c_str());
      mNCaloNormJets[itrig][iy] = mCaloDir.make<TH1F>(name,name,3,0,3);
      mNCaloNormJets[itrig][iy]->SetBit(TH1::kCanRebin);
      mNCaloNormJets[itrig][iy]->Sumw2();
      //-------- pt histograms ----------------------------
      if (mIsMC) {
        sprintf(name,"GenPt_Ybin%d",iy);
        mGenPt[iy] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxPt);
        mGenPt[iy]->Sumw2();
        sprintf(name,"GenX_Ybin%d",iy);
        mGenX[iy] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxX);
        mGenX[iy]->Sumw2();
      }
      sprintf(name,"PtPU_Ybin%d%s",iy,ss.c_str());
      mPFPtPU[itrig][iy] = mPFDir.make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mPFPtPU[itrig][iy]->Sumw2(); 
      sprintf(name,"Pt_Ybin%d%s",iy,ss.c_str());
      mPFPt[itrig][iy] = mPFDir.make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mPFPt[itrig][iy]->Sumw2();
      sprintf(name,"PtVsNPV_Ybin%d%s",iy,ss.c_str());
      mPFPtVsNPV[itrig][iy] = mPFDir.make<TH2F>(name,name,30,0,30,mPTBND.size()-1,auxPt);
      mPFPtVsNPV[itrig][iy]->Sumw2();
      sprintf(name,"NormPt_Ybin%d%s",iy,ss.c_str());
      mPFNormPt[itrig][iy] = mPFDir.make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mPFNormPt[itrig][iy]->Sumw2();
      sprintf(name,"PtVsNPV_Ybin%d%s",iy,ss.c_str());
      mCaloPtVsNPV[itrig][iy] = mCaloDir.make<TH2F>(name,name,30,0,30,mPTBND.size()-1,auxPt);
      mCaloPtVsNPV[itrig][iy]->Sumw2();
      sprintf(name,"Pt_Ybin%d%s",iy,ss.c_str());
      mCaloPt[itrig][iy] = mCaloDir.make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mCaloPt[itrig][iy]->Sumw2();
      sprintf(name,"NormPt_Ybin%d%s",iy,ss.c_str());
      mCaloNormPt[itrig][iy] = mCaloDir.make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mCaloNormPt[itrig][iy]->Sumw2();
      //-------- x histograms -----------------------------
      sprintf(name,"XEvents_Ybin%d%s",iy,ss.c_str());
      mPFEventsX[itrig][iy] = mPFDir.make<TH1F>(name,name,mPTBND.size()-1,auxX);
      mPFEventsX[itrig][iy]->Sumw2();
      sprintf(name,"Aux_Ybin%d%s",iy,ss.c_str());
      hAux[itrig][iy] = new TH1F(name,name,mPTBND.size()-1,auxX);
      sprintf(name,"PtVsNPV_Ybin%d%s",iy,ss.c_str());
      sprintf(name,"X_Ybin%d%s",iy,ss.c_str());
      mPFX[itrig][iy] = mPFDir.make<TH1F>(name,name,mPTBND.size()-1,auxX);
      mPFX[itrig][iy]->Sumw2();
      sprintf(name,"NormX_Ybin%d%s",iy,ss.c_str());
      mPFNormX[itrig][iy] = mPFDir.make<TH1F>(name,name,mPTBND.size()-1,auxX);
      mPFNormX[itrig][iy]->Sumw2();
      sprintf(name,"X_Ybin%d%s",iy,ss.c_str());
      mCaloX[itrig][iy] = mCaloDir.make<TH1F>(name,name,mPTBND.size()-1,auxX);
      mCaloX[itrig][iy]->Sumw2();
      sprintf(name,"NormX_Ybin%d%s",iy,ss.c_str());
      mCaloNormX[itrig][iy] = mCaloDir.make<TH1F>(name,name,mPTBND.size()-1,auxX);
      mCaloNormX[itrig][iy]->Sumw2();
      //-------- jet properties ----------------------------
      sprintf(name,"BetaStar_Ybin%d%s",iy,ss.c_str());
      mBetaStar[itrig][iy] = mPFDir.make<TH1F>(name,name,100,0,1.001);
      mBetaStar[itrig][iy]->Sumw2();
      sprintf(name,"CHF_Ybin%d%s",iy,ss.c_str());
      mCHF[itrig][iy] = mPFDir.make<TH1F>(name,name,100,0,1.001);
      mCHF[itrig][iy]->Sumw2();  
      sprintf(name,"NHF_Ybin%d%s",iy,ss.c_str());
      mNHF[itrig][iy] = mPFDir.make<TH1F>(name,name,100,0,1.001);
      mNHF[itrig][iy]->Sumw2();
      sprintf(name,"PHF_Ybin%d%s",iy,ss.c_str());
      mPHF[itrig][iy] = mPFDir.make<TH1F>(name,name,100,0,1.001);
      mPHF[itrig][iy]->Sumw2();
      sprintf(name,"ELF_Ybin%d%s",iy,ss.c_str());
      mELF[itrig][iy] = mPFDir.make<TH1F>(name,name,100,0,1.001);
      mELF[itrig][iy]->Sumw2();
      sprintf(name,"MUF_Ybin%d%s",iy,ss.c_str());
      mMUF[itrig][iy] = mPFDir.make<TH1F>(name,name,100,0,1.001);
      mMUF[itrig][iy]->Sumw2();
      sprintf(name,"JEC_Ybin%d%s",iy,ss.c_str());
      mJEC[itrig][iy] = mPFDir.make<TH1F>(name,name,500,0,10);
      mJEC[itrig][iy]->Sumw2();     

      sprintf(name,"CHFVsRun_Ybin%d%s",iy,ss.c_str()); 
      mCHFVsRun[itrig][iy] = mPFDir.make<TProfile>(name,name,3,0,3,0,1.1);
      mCHFVsRun[itrig][iy]->SetBit(TH1::kCanRebin);
      sprintf(name,"NHFVsRun_Ybin%d%s",iy,ss.c_str()); 
      mNHFVsRun[itrig][iy] = mPFDir.make<TProfile>(name,name,3,0,3,0,1.1);
      mNHFVsRun[itrig][iy]->SetBit(TH1::kCanRebin);
      sprintf(name,"PHFVsRun_Ybin%d%s",iy,ss.c_str()); 
      mPHFVsRun[itrig][iy] = mPFDir.make<TProfile>(name,name,3,0,3,0,1.1);
      mPHFVsRun[itrig][iy]->SetBit(TH1::kCanRebin);
      sprintf(name,"ELFVsRun_Ybin%d%s",iy,ss.c_str());
      mELFVsRun[itrig][iy] = mPFDir.make<TProfile>(name,name,3,0,3,0,1.1);
      mELFVsRun[itrig][iy]->SetBit(TH1::kCanRebin); 
      sprintf(name,"MUFVsRun_Ybin%d%s",iy,ss.c_str());
      mMUFVsRun[itrig][iy] = mPFDir.make<TProfile>(name,name,3,0,3,0,1.1);
      mMUFVsRun[itrig][iy]->SetBit(TH1::kCanRebin);
     
      sprintf(name,"EMF_Ybin%d%s",iy,ss.c_str());
      mEMF[itrig][iy] = mCaloDir.make<TH1F>(name,name,100,0,1.001);
      mEMF[itrig][iy]->Sumw2();
      sprintf(name,"N90hits_Ybin%d%s",iy,ss.c_str());
      mN90hits[itrig][iy] = mCaloDir.make<TH1F>(name,name,200,0,200);
      mN90hits[itrig][iy]->Sumw2();
      sprintf(name,"fHPD_Ybin%d%s",iy,ss.c_str());
      mfHPD[itrig][iy] = mCaloDir.make<TH1F>(name,name,100,0,1.001);
      mfHPD[itrig][iy]->Sumw2();
      sprintf(name,"NTrkCalo_Ybin%d%s",iy,ss.c_str());
      mNTrkCalo[itrig][iy] = mCaloDir.make<TH1F>(name,name,200,0,200);
      mNTrkCalo[itrig][iy]->Sumw2();
      sprintf(name,"NTrkVtx_Ybin%d%s",iy,ss.c_str());
      mNTrkVtx[itrig][iy] = mCaloDir.make<TH1F>(name,name,200,0,200);
      mNTrkVtx[itrig][iy]->Sumw2();

      sprintf(name,"EMFVsRun_Ybin%d%s",iy,ss.c_str()); 
      mEMFVsRun[itrig][iy] = mCaloDir.make<TProfile>(name,name,3,0,3,0,1.1);
      mEMFVsRun[itrig][iy]->SetBit(TH1::kCanRebin);
      sprintf(name,"NTrkCaloVsRun_Ybin%d%s",iy,ss.c_str());
      mNTrkCaloVsRun[itrig][iy] = mCaloDir.make<TProfile>(name,name,3,0,3,0,1.1);
      mNTrkCaloVsRun[itrig][iy]->SetBit(TH1::kCanRebin);
      sprintf(name,"NTrkVtxVsRun_Ybin%d%s",iy,ss.c_str());
      mNTrkVtxVsRun[itrig][iy] = mCaloDir.make<TProfile>(name,name,3,0,3,0,1.1);
      mNTrkVtxVsRun[itrig][iy]->SetBit(TH1::kCanRebin);
    }// y loop
  }// trigger loop
}
//////////////////////////////////////////////////////////////////////////////////////////
void InclusiveHistos::endJob() 
{

}
//////////////////////////////////////////////////////////////////////////////////////////
int InclusiveHistos::getBin(double x, const std::vector<double>& boundaries)
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
int InclusiveHistos::findRun(int x, const std::vector<int>& runs)
{
  int result(-1);
  for(unsigned i=0;i<runs.size();i++)
    if (x == runs[i])
      return i;
  return result;
}
//////////////////////////////////////////////////////////////////////////////////////////
void InclusiveHistos::analyze(edm::Event const& evt, edm::EventSetup const& iSetup) 
{ 
  ofstream logFile;
  logFile.open(mLogName.c_str());
  logFile.setf(ios::right);
  int decade = 0;
  int counter_hlt[100],counter_pv(0),counter_hcal(0);
  int pf_counter_y[100][10],pf_counter_pt[100][10],pf_counter_id[100][10],pf_counter_met[100],pf_counter_rho[100];
  int calo_counter_y[100][10],calo_counter_pt[100][10],calo_counter_id[100][10],calo_counter_met[100],calo_counter_rho[100];
  int Ntrig = (int)mTriggers.size();
  vector<int> vRuns(0);
  for(int itrig=0;itrig<TMath::Max(1,Ntrig);itrig++) {
    counter_hlt[itrig] = 0;
    pf_counter_met[itrig] = 0;
    pf_counter_rho[itrig] = 0;
    calo_counter_met[itrig] = 0;
    calo_counter_rho[itrig] = 0;
    for(int iy=0;iy<10;iy++) {
      pf_counter_y[itrig][iy]    = 0;
      pf_counter_pt[itrig][iy]   = 0;
      pf_counter_id[itrig][iy]   = 0;
      calo_counter_y[itrig][iy]  = 0;
      calo_counter_pt[itrig][iy] = 0;
      calo_counter_id[itrig][iy] = 0;
    }
  }
  //---------- loop over the files -----------------------
  for(unsigned iFile=0;iFile<mFileNames.size();iFile++) {
    TFile *inf = TFile::Open(mFileNames[iFile].c_str());  
    TTree *mTree = (TTree*)inf->Get((mDirName+"/"+mTreeName).c_str());
    unsigned long NEntries = mTree->GetEntries();
    cout<<"Adding file: "<<mFileNames[iFile]<<endl;
    cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
    logFile<<"Adding file: "<<mFileNames[iFile]<<"\n";
    logFile<<"Reading TREE: "<<NEntries<<" events"<<"\n";   
    TBranch *branch = mTree->GetBranch("events");
    branch->SetAddress(&mEvent);
    //---------- loop over the events ----------------------
    unsigned long NN = NEntries;
    if (mNEvents > -1)
      NN = (unsigned long)mNEvents;
    cout<<"Will run over "<<NN<<" events"<<endl;
    for(unsigned long i=0;i<NN;i++) {
      double progress = 10.0*i/(1.0*NN);
      int k = TMath::FloorNint(progress); 
      if (k > decade) 
        cout<<10*k<<" %"<<endl;
      decade = k;          
      mTree->GetEntry(i);
      double wt(1.0);
      char name[1000],runName[1000];
      if (mIsMC) {
        wt = mEvent->evtHdr().weight();
        mGenPU->Fill(mEvent->evtHdr().intpu());
        if (mPUFileName != "" && mPUHistName != "") {
          int nINTpu = mEvent->evtHdr().intpu();
          int nOOTpu = mEvent->evtHdr().ootpuEarly()+mEvent->evtHdr().ootpuLate();
          double wtINT(0.0),wtOOT(0.0);
          if (nINTpu < mPUh->GetNbinsX())   
            wtINT = mPUh->GetBinContent(nINTpu+1); 
          if (nOOTpu < mPUh->GetNbinsX())
            wtOOT = mPUh->GetBinContent(nOOTpu+1);
          wt *= wtINT;
          //cout<<"INTpu = "<<nINTpu<<", wINT = "<<wtINT<<", OOTpu = "<<nOOTpu<<", wOOT = "<<wtOOT<<", total weight = "<<wtINT*wtOOT<<endl;
        }
        for(unsigned j=0;j<mEvent->nGenJets();j++) {
          int ybin = getBin(fabs(mEvent->genjet(j).Rapidity()),mYBND);
          if (ybin > -1) {
            mGenPt[ybin]->Fill(mEvent->genjet(j).pt(),wt);
            mGenX[ybin]->Fill(0.001*mEvent->genjet(j).pt(),wt);
          }
        }
      }
      sprintf(runName,"%d",mEvent->evtHdr().runNo());
      mPFRhoVsRun->Fill(runName,mEvent->evtHdr().pfRho(),wt);
      mPFRhoVsNPV->Fill(mEvent->evtHdr().nVtxGood(),mEvent->evtHdr().pfRho(),wt);
      mCaloRhoVsRun->Fill(runName,mEvent->evtHdr().caloRho(),wt);
      mCaloRhoVsNPV->Fill(mEvent->evtHdr().nVtxGood(),mEvent->evtHdr().caloRho(),wt);
      mNPV->Fill(mEvent->evtHdr().nVtx(),wt);
      mPVx->Fill(mEvent->evtHdr().PVx(),wt);
      mPVy->Fill(mEvent->evtHdr().PVy(),wt);
      mPVz->Fill(mEvent->evtHdr().PVz(),wt);
      mBSx->Fill(mEvent->evtHdr().BSx(),wt);
      mBSy->Fill(mEvent->evtHdr().BSy(),wt);
      mBSz->Fill(mEvent->evtHdr().BSz(),wt);
      if (!mIsMC) {
        mNvtxVsRun->Fill(runName,mEvent->evtHdr().nVtxGood(),wt);
      }
      if (!mEvent->evtHdr().isPVgood()) continue; 
      counter_pv++;
      bool HBHEfilter(true);
      if (mApplyHBEHfilter) {
        HBHEfilter = (mEvent->evtHdr().hcalNoise() == 1);
      } 
      if (!HBHEfilter) continue;
      counter_hcal++;
      //---------- loop over the triggers ----------------------
      for(int itrig=0;itrig<TMath::Max(1,Ntrig);itrig++) {
        int prescale(1),ihlt;
        bool hltPass(false);
        if (Ntrig == 0) {
          ihlt = 0;
          hltPass = true;
        }
        else {
          //-------- loop over the versions ----------------------
          unsigned nVersions = mTrigIndex[itrig].size();
          for(unsigned iv=0;iv<nVersions;iv++) {
            ihlt = mTrigIndex[itrig][iv];
            if (mEvent->fired(ihlt) > 0) {
              hltPass = true;
              prescale = mEvent->preL1(ihlt) * mEvent->preHLT(ihlt);
              continue;
            }
          }
        }
        //---- reset the auxilary file for event counting ------
        for(unsigned iy=0;iy<mYBND.size()-1;iy++) {
          hAux[itrig][iy]->Reset();
        }
        //cout<<ihlt<<" "<<mEvent->fired(ihlt)<<" "<<hltPass<<endl;
        //-------- cut flow ------------------------------------
        if (hltPass) {
          counter_hlt[itrig]++;
          int nPFGoodJets(0);
          int nCaloGoodJets(0);
          bool cut_PFMET = (mEvent->pfmet().met_o_sumet() < mMaxMETovSumET);
          if (cut_PFMET) {
            pf_counter_met[itrig]++;
            bool cut_PFRHO = (mEvent->evtHdr().pfRho() < 100);
            if (cut_PFRHO) {
              pf_counter_rho[itrig]++;
              for(unsigned j=0;j<mEvent->nPFJets();j++) {
                double y = mEvent->pfjet(j).y();
                int ybin = getBin(fabs(y),mYBND);  
                if (ybin < 0) continue;
                pf_counter_y[itrig][ybin]++;
                double ptCor = mEvent->pfjet(j).ptCor();
                bool cut_pt = (mEvent->pfjet(j).ptCor() >= mMinPt[itrig]);
                if (!cut_pt) continue;
                pf_counter_pt[itrig][ybin]++;
                bool cutID=(mEvent->pfjet(j).tightID() && (mEvent->pfjet(j).elf() < 0.9) && (mEvent->pfjet(j).muf() < 0.9) && (mEvent->pfjet(j).nhf() < 0.9));
                if (!cutID) continue;
                bool cut_BetaStar = ((fabs(y) > 2.5) || (mEvent->pfjet(j).betaStar() < mMaxBetaStar));
                if (!cut_BetaStar) {
                  mPFPtPU[itrig][ybin]->Fill(ptCor,wt);
                  continue;
                }
                pf_counter_id[itrig][ybin]++;
                nPFGoodJets++;
                mPFPt[itrig][ybin]->Fill(ptCor,wt);
                mPFPtVsNPV[itrig][ybin]->Fill(mEvent->evtHdr().nVtxGood(),ptCor,prescale);
                mPFNormPt[itrig][ybin]->Fill(ptCor,prescale);
                double x = ptCor/1000.;
                // ----- fill the independent events ----------
                int xbin = hAux[itrig][ybin]->FindBin(x);
                double z = hAux[itrig][ybin]->GetBinContent(xbin);
                if (z == 0) {
                  hAux[itrig][ybin]->Fill(x);
                }
                mPFX[itrig][ybin]->Fill(x,wt);
                mPFNormX[itrig][ybin]->Fill(x,prescale);
                mBetaStar[itrig][ybin]->Fill(mEvent->pfjet(j).betaStar(),wt);
                mCHF[itrig][ybin]->Fill(mEvent->pfjet(j).chf(),wt);
                mNHF[itrig][ybin]->Fill(mEvent->pfjet(j).nhf(),wt);
                mPHF[itrig][ybin]->Fill(mEvent->pfjet(j).phf(),wt); 
                mELF[itrig][ybin]->Fill(mEvent->pfjet(j).elf(),wt);
                mMUF[itrig][ybin]->Fill(mEvent->pfjet(j).muf(),wt);
                mJEC[itrig][ybin]->Fill(mEvent->pfjet(j).cor(),wt);
                mNPFJets[itrig][ybin]->Fill(runName,wt);
                mNPFNormJets[itrig][ybin]->Fill(runName,prescale);
                mCHFVsRun[itrig][ybin]->Fill(runName,mEvent->pfjet(j).chf(),wt);
                mNHFVsRun[itrig][ybin]->Fill(runName,mEvent->pfjet(j).nhf(),wt); 
                mPHFVsRun[itrig][ybin]->Fill(runName,mEvent->pfjet(j).phf(),wt);
                mELFVsRun[itrig][ybin]->Fill(runName,mEvent->pfjet(j).elf(),wt);
                mMUFVsRun[itrig][ybin]->Fill(runName,mEvent->pfjet(j).muf(),wt);
              }// pfjet loop  
              mPFJetMulti[itrig]->Fill(nPFGoodJets,wt); 
              if (nPFGoodJets > 0) {
                mPFMETovSUMET[itrig]->Fill(mEvent->pfmet().met_o_sumet(),wt);
              }
            }// pf rho cut
          }// pf met cut
          for(unsigned iy=0;iy<mYBND.size()-1;iy++) {
            mPFEventsX[itrig][iy]->Add(hAux[itrig][iy]); 
          }
          bool cut_CaloMET = (mEvent->calomet().met_o_sumet() < mMaxMETovSumET);
          if (cut_CaloMET) {
            calo_counter_met[itrig]++;
            bool cut_CALORHO = (mEvent->evtHdr().caloRho() < 100);
            if (cut_CALORHO) {
              calo_counter_rho[itrig]++;
              for(unsigned j=0;j<mEvent->nCaloJets();j++) {
                double y = mEvent->calojet(j).y();
                int ybin = getBin(fabs(y),mYBND);
                if (ybin < 0) continue;
                calo_counter_y[itrig][ybin]++;
                double ptCor = mEvent->calojet(j).ptCor();
                bool cut_pt = (ptCor >= mMinPt[itrig]);
                if (!cut_pt) continue;
                calo_counter_pt[itrig][ybin]++;
                bool cutID = (mEvent->calojet(j).tightID() && (mEvent->calojet(j).emf() > 0.05));
                if (!cutID) continue;
                calo_counter_id[itrig][ybin]++;
                nCaloGoodJets++;
                mCaloPt[itrig][ybin]->Fill(ptCor,wt);
                mCaloPtVsNPV[itrig][ybin]->Fill(mEvent->evtHdr().nVtxGood(),ptCor,prescale);
                mCaloNormPt[itrig][ybin]->Fill(ptCor,prescale);
                double x = ptCor/1000.;
                mCaloX[itrig][ybin]->Fill(x,wt);
                mCaloNormX[itrig][ybin]->Fill(x,prescale);
                mEMF[itrig][ybin] ->Fill((mEvent->calojet(j)).emf(),wt);
                mN90hits[itrig][ybin]->Fill((mEvent->calojet(j)).n90hits(),wt);
                mfHPD[itrig][ybin]->Fill((mEvent->calojet(j)).fHPD(),wt);
                mNTrkCalo[itrig][ybin]->Fill((mEvent->calojet(j)).nTrkCalo(),wt);
                mNTrkVtx[itrig][ybin]->Fill((mEvent->calojet(j)).nTrkVtx(),wt);              
                mNCaloJets[itrig][ybin]->Fill(runName,wt);
                mNCaloNormJets[itrig][ybin]->Fill(runName,prescale);
                mEMFVsRun[itrig][ybin]->Fill(name,(mEvent->calojet(j)).emf(),wt);
                mNTrkCaloVsRun[itrig][ybin]->Fill(runName,(mEvent->calojet(j)).nTrkCalo(),wt);
                mNTrkVtxVsRun[itrig][ybin]->Fill(runName,(mEvent->calojet(j)).nTrkVtx(),wt);
              }// calojet loop
              mCaloJetMulti[itrig]->Fill(nCaloGoodJets,wt);
              if (nCaloGoodJets > 0) {
                mCaloMETovSUMET[itrig]->Fill(mEvent->calomet().met_o_sumet(),wt);
              }
            } // calo rho cut
          } // calo met cut  
        }// if hlt
      }// trigger loop
    }// tree loop
  }// file loop
  logFile<<"Events passing the PV:         "<<counter_pv<<"\n";
  logFile<<"Events passing the HBHE noise:         "<<counter_hcal<<"\n";
  for(int itrig=0;itrig<TMath::Max(1,Ntrig);itrig++) {
    logFile<<"*********************************************"<<"\n";
    if (Ntrig > 0) {
      logFile<<"Trigger path: "              <<mTriggers[itrig]<<"\n";
      logFile<<"Events passing the trigger:    "<<counter_hlt[itrig]<<"\n";
    } 
    logFile<<"Events passing the PFMET:      "<<pf_counter_met[itrig]<<"\n";
    logFile<<"Events passing the PFRHO:      "<<pf_counter_rho[itrig]<<"\n";
    logFile<<"Events passing the CaloMET:    "<<calo_counter_met[itrig]<<"\n";
    logFile<<"Events passing the CaloRHO:    "<<calo_counter_rho[itrig]<<"\n";
    for(unsigned j=0;j<mYBND.size()-1;j++) {
      logFile<<"--------------------------------------------"<<"\n"; 
      logFile<<"["<<mYBND[j]<<", "<<mYBND[j+1]<<"]"<<"\n";
      logFile<<"pt > "<<mMinPt[itrig]<<" GeV, PFJets =  "<<pf_counter_pt[itrig][j]<<", CaloJets = "<<calo_counter_pt[itrig][j]<<"\n";
      logFile<<"id, PFJets = "<<pf_counter_id[itrig][j]<<", CaloJets = "<<calo_counter_id[itrig][j]<<"\n";
    }
  }
  logFile.close();
}
//////////////////////////////////////////////////////////////////////////////////////////
InclusiveHistos::~InclusiveHistos() 
{
}

DEFINE_FWK_MODULE(InclusiveHistos);
