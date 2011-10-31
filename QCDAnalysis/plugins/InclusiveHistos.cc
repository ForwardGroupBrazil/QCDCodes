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
  mYBND       = cfg.getParameter<std::vector<double> > ("yBnd");
  mPTBND      = cfg.getParameter<std::vector<double> > ("ptBnd");
  mMinPt      = cfg.getParameter<std::vector<double> > ("minPt");
  mFileNames  = cfg.getParameter<std::vector<std::string> > ("filenames");
  mTreeName   = cfg.getParameter<std::string> ("treename");
  mDirName    = cfg.getParameter<std::string> ("dirname");
  mLogName    = cfg.getParameter<std::string> ("logname");
  mPUFileName = cfg.getUntrackedParameter<std::string> ("puFile","");
  mPUHistName = cfg.getUntrackedParameter<std::string> ("puHisto","");
  mTriggers   = cfg.getParameter<std::vector<std::string> > ("triggers");
  mIsMC       = cfg.getParameter<bool> ("isMC");
  mNEvents    = cfg.getParameter<int>  ("nEvents"); 
  mJetID      = cfg.getParameter<int>  ("jetID");
  mHCALNoise  = cfg.getParameter<int>  ("hcalNoiseFilter");
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
  /*
  mTree = new TChain((mDirName+"/"+mTreeName).c_str());
  for(unsigned ifile=0;ifile < mFileNames.size();ifile++) {
    cout<<"Adding file: "<<mFileNames[ifile]<<endl;
    mTree->Add((mFileNames[ifile]).c_str(),-1);
    cout<<"entries: "<<mTree->GetEntries()<<endl;
  }
  TBranch *branch = mTree->GetBranch("events");
  branch->SetAddress(&mEvent);
  */
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
    for(unsigned itrig=0;itrig<mTriggers.size();itrig++) {
      int index(-1);
      for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
        string ss = hTrigNames->GetXaxis()->GetBinLabel(ibin+1);
        if (ss == mTriggers[itrig]) {
          index = ibin;
          continue;
        }
      }
      if (index > -1)
        mTrigIndex.push_back(index); 
      else {
        throw cms::Exception("InclusiveHistos: ")<<"The requested trigger ("<<mTriggers[itrig]<<") is not found ";
      }
    }
    for(unsigned itrig=0;itrig<mTriggers.size();itrig++)
      cout<<mTriggers[itrig]<<" --> "<<mTrigIndex[itrig]<<endl;
  } 
  //--------- book histos ---------------------------------
  TFileDirectory mPFDir   = fs->mkdir(mDirName+"pf");
  TFileDirectory mCaloDir = fs->mkdir(mDirName+"calo");
  double auxPt[500],auxX[500];
  for(unsigned ipt=0;ipt<mPTBND.size();ipt++) {
    auxPt[ipt] = mPTBND[ipt];
    auxX[ipt] = mPTBND[ipt]/1000.; 
  }
  mNPV = fs->make<TH1F>("NPV","NPV",30,0,30);
  mPVx = fs->make<TH1F>("PVx","PVx",200,-0.1,0.1);
  mPVy = fs->make<TH1F>("PVy","PVy",200,-0.1,0.1);
  mPVz = fs->make<TH1F>("PVz","PVz",200,-30,30);
  mBSx = fs->make<TH1F>("BSx","BSx",200,-0.1,0.1);
  mBSy = fs->make<TH1F>("BSy","BSy",200,-0.1,0.1);
  mBSz = fs->make<TH1F>("BSz","BSz",200,-0.2,0.2);
  mCaloRhoVsRun = fs->make<TProfile>("CaloRhoVsRun","CaloRhoVsRun",3,0,3,0,50);
  mCaloRhoVsRun->SetBit(TH1::kCanRebin);
  mPFRhoVsRun = fs->make<TProfile>("PFRhoVsRun","PFRhoVsRun",3,0,3,0,50);
  mPFRhoVsRun->SetBit(TH1::kCanRebin);
  mCaloRhoVsNPV = fs->make<TProfile>("CaloRhoVsNPV","CaloRhoVsNPV",30,0,30,0,50);
  mPFRhoVsNPV = fs->make<TProfile>("PFRhoVsNPV","PFRhoVsNPV",30,0,30,0,50);
  int Ntrig = (int)mTriggers.size();
  for(int itrig=0;itrig<TMath::Max(Ntrig,1);itrig++) {
    string ss("");
    if (Ntrig > 0)
      ss = "_"+mTriggers[itrig];
    for(int im=0;im<4;im++) {
      sprintf(name,"METovSUMET_%dJet%s",im+1,ss.c_str());
      mPFMETovSUMET[itrig][im]   = mPFDir.make<TH1F>(name,name,100,0,1.0001);
      mCaloMETovSUMET[itrig][im] = mCaloDir.make<TH1F>(name,name,100,0,1.0001);
    }
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
  int counter_hlt[100],counter_pv[100],counter_hcal[100];
  int pf_counter_y[100][10],pf_counter_pt[100][10],pf_counter_id[100][10];
  int calo_counter_y[100][10],calo_counter_pt[100][10],calo_counter_id[100][10];
  int Ntrig = (int)mTriggers.size();
  vector<int> Runs;
  for(int itrig=0;itrig<TMath::Max(1,Ntrig);itrig++) {
    counter_hlt[itrig]  = 0;
    counter_pv[itrig]   = 0;
    counter_hcal[itrig] = 0;
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
      char name[1000];
      if (mIsMC) {
        //wt = mEvent->evtHdr().weight();
        if (mPUFileName != "" && mPUHistName != "") {
          int nINTpu = mEvent->evtHdr().intpu();
          int nOOTpu = mEvent->evtHdr().ootpuEarly()+mEvent->evtHdr().ootpuLate();
          double wtINT(0.0),wtOOT(0.0);
          if (nINTpu < mPUh->GetNbinsX())   
            wtINT = mPUh->GetBinContent(nINTpu+1); 
          if (nOOTpu < mPUh->GetNbinsX())
            wtOOT = mPUh->GetBinContent(nOOTpu+1);
          wt *= wtINT*wtOOT;
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
      //---------- loop over the triggers ----------------------
      for(int itrig=0;itrig<TMath::Max(1,Ntrig);itrig++) {
        int prescale(1),ihlt;
        bool hltPass(false);
        if (Ntrig == 0) {
          ihlt = 0;
          hltPass = true;
        }
        else {
          ihlt = mTrigIndex[itrig];
          if (mEvent->fired(ihlt) > 0) {
            hltPass = true;
            prescale = mEvent->preL1(ihlt) * mEvent->preHLT(ihlt);
          }
        }
        //cout<<ihlt<<" "<<mEvent->fired(ihlt)<<" "<<hltPass<<endl;
        //-------- cut flow ------------------------------------
        if (hltPass) {
          counter_hlt[itrig]++;
          if (mEvent->evtHdr().isPVgood() == 1) {
            counter_pv[itrig]++;
            bool cut_hcalNoise(true);
            if (mHCALNoise == 1)
              cut_hcalNoise = mEvent->evtHdr().looseHCALNoise();
            if (mHCALNoise == 2)
              cut_hcalNoise = mEvent->evtHdr().tightHCALNoise(); 
            if (cut_hcalNoise) {
              counter_hcal[itrig]++;  
              int nPFGoodJets(0);
              for(unsigned j=0;j<mEvent->nPFJets();j++) {
                double y = mEvent->pfjet(j).y();
                int ybin = getBin(fabs(y),mYBND);  
                if (ybin > -1) {
                  pf_counter_y[itrig][ybin]++;
                  if (mEvent->pfjet(j).ptCor() >= mMinPt[itrig]) {
                    pf_counter_pt[itrig][ybin]++;
                    bool cutID(true);
                    if (mJetID == 1)
                      cutID = mEvent->pfjet(j).looseID();
                    if (mJetID == 2)
                      cutID = mEvent->pfjet(j).tightID();  
                    if (cutID) {
                      pf_counter_id[itrig][ybin]++;
                      nPFGoodJets++;
                      double ptCor = mEvent->pfjet(j).ptCor();
                      mPFPt[itrig][ybin]->Fill(ptCor,wt);
                      mPFPtVsNPV[itrig][ybin]->Fill(mEvent->evtHdr().nVtxGood(),ptCor,prescale);
                      mPFNormPt[itrig][ybin]->Fill(ptCor,prescale);
                      double x = ptCor/1000.;
                      mPFX[itrig][ybin]->Fill(x,wt);
                      mPFNormX[itrig][ybin]->Fill(x,prescale);
                      mCHF[itrig][ybin]->Fill((mEvent->pfjet(j)).chf(),wt);
                      mNHF[itrig][ybin]->Fill((mEvent->pfjet(j)).nhf(),wt);
                      mPHF[itrig][ybin]->Fill((mEvent->pfjet(j)).phf(),wt); 
                      mELF[itrig][ybin]->Fill((mEvent->pfjet(j)).elf(),wt);
                      sprintf(name,"%d",mEvent->evtHdr().runNo());	
                      mNPFJets[itrig][ybin]->Fill(name,1);
                      mNPFNormJets[itrig][ybin]->Fill(name,prescale);
                      mCHFVsRun[itrig][ybin]->Fill(name,(mEvent->pfjet(j)).chf(),1);
                      mNHFVsRun[itrig][ybin]->Fill(name,(mEvent->pfjet(j)).nhf(),1); 
                      mPHFVsRun[itrig][ybin]->Fill(name,(mEvent->pfjet(j)).phf(),1);
                      mELFVsRun[itrig][ybin]->Fill(name,(mEvent->pfjet(j)).elf(),1);
                    }// cut id
                  }// cut max pt
                }// cut y  
              }// pfjet loop
              if (nPFGoodJets > 0) {
                sprintf(name,"%d",mEvent->evtHdr().runNo());
                mPFRhoVsRun->Fill(name,mEvent->evtHdr().pfRho(),1);
                mPFRhoVsNPV->Fill(mEvent->evtHdr().nVtxGood(),mEvent->evtHdr().pfRho(),1);
              }
              mPFJetMulti[itrig]->Fill(nPFGoodJets,wt);
              for(int ii=0;ii<4;ii++) {
                if (nPFGoodJets == ii+1)
                  mPFMETovSUMET[itrig][ii]->Fill(mEvent->pfmet().met_o_sumet(),wt);
              }
              int nCaloGoodJets(0);
              for(unsigned j=0;j<mEvent->nCaloJets();j++) {
                double y = mEvent->calojet(j).y();
                int ybin = getBin(fabs(y),mYBND);
                if (ybin > -1) {
                  calo_counter_y[itrig][ybin]++;
                  if (mEvent->calojet(j).ptCor() >= mMinPt[itrig]) {
                    calo_counter_pt[itrig][ybin]++;
                    bool cutID(true);
                    if (mJetID == 1)
                      cutID = mEvent->calojet(j).looseID();
                    if (mJetID == 2)
                      cutID = mEvent->calojet(j).tightID();
                    if (cutID) {
                      calo_counter_id[itrig][ybin]++;
                      nCaloGoodJets++;
                      double ptCor = mEvent->calojet(j).ptCor();
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
                      sprintf(name,"%d",mEvent->evtHdr().runNo());
                      mNCaloJets[itrig][ybin]->Fill(name,1);
                      mNCaloNormJets[itrig][ybin]->Fill(name,prescale);
                      mEMFVsRun[itrig][ybin]->Fill(name,(mEvent->calojet(j)).emf(),1);
                      mNTrkCaloVsRun[itrig][ybin]->Fill(name,(mEvent->calojet(j)).nTrkCalo(),1);
                      mNTrkVtxVsRun[itrig][ybin]->Fill(name,(mEvent->calojet(j)).nTrkVtx(),1);
                    }// cut id
                  }// cut min pt
                }// cut y
              }// calojet loop
              mCaloJetMulti[itrig]->Fill(nCaloGoodJets,wt);
              for(int ii=0;ii<4;ii++) {
                if (nCaloGoodJets == ii+1)
                  mCaloMETovSUMET[itrig][ii]->Fill(mEvent->calomet().met_o_sumet(),wt);
              }
              if (nCaloGoodJets > 0) {
                sprintf(name,"%d",mEvent->evtHdr().runNo());
                mCaloRhoVsRun->Fill(name,mEvent->evtHdr().caloRho(),1);
                mCaloRhoVsNPV->Fill(mEvent->evtHdr().nVtxGood(),mEvent->evtHdr().caloRho(),1);
              } 
              if (nPFGoodJets > 0 && nCaloGoodJets > 0) {
                mNPV->Fill(mEvent->evtHdr().nVtx());
                mPVx->Fill(mEvent->evtHdr().PVx());
                mPVy->Fill(mEvent->evtHdr().PVy());
                mPVz->Fill(mEvent->evtHdr().PVz());
                mBSx->Fill(mEvent->evtHdr().BSx());
                mBSy->Fill(mEvent->evtHdr().BSy());
                mBSz->Fill(mEvent->evtHdr().BSz());
              }
            }// if hcal noise
          }// if pv
        }// if hlt
      }// trigger loop
    }// tree loop
  }// file loop
  for(int itrig=0;itrig<TMath::Max(1,Ntrig);itrig++) {
    logFile<<"*********************************************"<<"\n";
    if (Ntrig > 0) {
      logFile<<"Trigger path: "              <<mTriggers[itrig]<<"\n";
      logFile<<"Events passing the trigger:    "<<counter_hlt[itrig]<<"\n";
    } 
    logFile<<"Events passing the PV:         "<<counter_pv[itrig]<<"\n";
    logFile<<"Events passing the HCAL NOISE: "<<counter_hcal[itrig]<<"\n";  
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
