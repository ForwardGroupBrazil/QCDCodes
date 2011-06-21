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
  mYBND      = cfg.getParameter<std::vector<double> > ("yBnd");
  mPTBND     = cfg.getParameter<std::vector<double> > ("ptBnd");
  mMinPt     = cfg.getParameter<std::vector<double> > ("minPt");
  mFileName  = cfg.getParameter<std::string> ("filename");
  mTreeName  = cfg.getParameter<std::string> ("treename");
  mDirName   = cfg.getParameter<std::string> ("dirname");
  mTriggers  = cfg.getParameter<std::vector<std::string> > ("triggers");
  mIsMC      = cfg.getParameter<bool> ("isMC");
  mNEvents   = cfg.getParameter<int>  ("nEvents"); 
  mJetID     = cfg.getParameter<int>  ("jetID");
  mHCALNoise = cfg.getParameter<int>  ("hcalNoiseFilter");
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
  mInf = TFile::Open(mFileName.c_str());
  mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());
  mTree = (TTree*)mDir->Get(mTreeName.c_str());
  mEvent = new QCDEvent();
  TBranch *branch = mTree->GetBranch("events");
  branch->SetAddress(&mEvent);
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
      sprintf(name,"NCaloJetsVsRun_Ybin%d%s",iy,ss.c_str());
      mNCaloJets[itrig][iy] = mCaloDir.make<TH1F>(name,name,3,0,3);
      mNCaloJets[itrig][iy]->SetBit(TH1::kCanRebin);
      mNCaloJets[itrig][iy]->Sumw2();
      //-------- pt histograms ----------------------------
      sprintf(name,"Pt_Ybin%d%s",iy,ss.c_str());
      mPFPt[itrig][iy] = mPFDir.make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mPFPt[itrig][iy]->Sumw2();
      sprintf(name,"NormPt_Ybin%d%s",iy,ss.c_str());
      mPFNormPt[itrig][iy] = mPFDir.make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mPFNormPt[itrig][iy]->Sumw2();
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
    }// y loop
  }// trigger loop
}
//////////////////////////////////////////////////////////////////////////////////////////
void InclusiveHistos::endJob() 
{
  mInf->Close();
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
  unsigned NEntries = mTree->GetEntries();
  cout<<"File: "<<mFileName<<endl;
  cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
  int decade = 0;
  int counter_hlt[50],counter_pv[50],counter_hcal[50];
  int pf_counter_y[50][10],pf_counter_pt[50][10],pf_counter_id[50][10];
  int calo_counter_y[50][10],calo_counter_pt[50][10],calo_counter_id[50][10];
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
  //---------- loop over the events ----------------------
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
    double wt(1.0);
    char name[1000];
    if (mIsMC)
      wt = mEvent->evtHdr().weight(); 
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
              int ybin = getBin(fabs(mEvent->pfjet(j).y()),mYBND);  
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
                    mPFPt[itrig][ybin]->Fill((mEvent->pfjet(j)).ptCor(),wt);
                    mPFNormPt[itrig][ybin]->Fill((mEvent->pfjet(j)).ptCor(),prescale);
                    double x = mEvent->pfjet(j).ptCor()/1000.;
                    mPFX[itrig][ybin]->Fill(x,wt);
                    mPFNormX[itrig][ybin]->Fill(x,prescale);
                    mCHF[itrig][ybin]->Fill((mEvent->pfjet(j)).chf(),wt);
                    mNHF[itrig][ybin]->Fill((mEvent->pfjet(j)).nhf(),wt);
                    mPHF[itrig][ybin]->Fill((mEvent->pfjet(j)).phf(),wt); 
                    sprintf(name,"%d",mEvent->evtHdr().runNo());	
                    mNPFJets[itrig][ybin]->Fill(name,prescale);
                  }// cut id
                }// cut max pt
              }// cut y  
            }// pfjet loop
            mPFJetMulti[itrig]->Fill(nPFGoodJets,wt);
            for(int ii=0;ii<4;ii++) {
              if (nPFGoodJets > ii)
                mPFMETovSUMET[itrig][ii]->Fill(mEvent->pfmet().met_o_sumet(),wt);
            }
            int nCaloGoodJets(0);
            for(unsigned j=0;j<mEvent->nCaloJets();j++) {
              int ybin = getBin(fabs(mEvent->calojet(j).y()),mYBND);
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
                    mCaloPt[itrig][ybin]->Fill((mEvent->calojet(j)).ptCor(),wt);
                    mCaloNormPt[itrig][ybin]->Fill((mEvent->calojet(j)).ptCor(),prescale);
                    double x = mEvent->calojet(j).ptCor()/1000.;
                    mCaloX[itrig][ybin]->Fill(x,wt);
                    mCaloNormX[itrig][ybin]->Fill(x,prescale);
                    mEMF[itrig][ybin] ->Fill((mEvent->calojet(j)).emf(),wt);
                    mN90hits[itrig][ybin]->Fill((mEvent->calojet(j)).n90hits(),wt);
                    mfHPD[itrig][ybin]->Fill((mEvent->calojet(j)).fHPD(),wt);
                    mNTrkCalo[itrig][ybin]->Fill((mEvent->calojet(j)).nTrkCalo(),wt);
                    mNTrkVtx[itrig][ybin]->Fill((mEvent->calojet(j)).nTrkVtx(),wt);              
                    sprintf(name,"%d",mEvent->evtHdr().runNo());
                    mNCaloJets[itrig][ybin]->Fill(name,prescale);
                  }// cut id
                }// cut min pt
              }// cut y
            }// calojet loop
            mCaloJetMulti[itrig]->Fill(nCaloGoodJets,wt);
            for(int ii=0;ii<4;ii++) {
              if (nCaloGoodJets > ii)
                mCaloMETovSUMET[itrig][ii]->Fill(mEvent->calomet().met_o_sumet(),wt);
            }
          }// if hcal noise
        }// if pv
      }// if hlt
    }// trigger loop
  }// tree loop
  for(int itrig=0;itrig<TMath::Max(1,Ntrig);itrig++) {
    cout<<"*********************************************"<<endl;
    if (Ntrig > 0) {
      cout<<"Trigger path: "              <<mTriggers[itrig]<<endl;
      cout<<"Events passing the trigger:    "<<counter_hlt[itrig]<<endl;
    } 
    cout<<"Events passing the PV:         "<<counter_pv[itrig]<<endl;
    cout<<"Events passing the HCAL NOISE: "<<counter_hcal[itrig]<<endl;  
    for(unsigned j=0;j<mYBND.size()-1;j++) {
      cout<<"--------------------------------------------"<<endl; 
      cout<<"["<<mYBND[j]<<", "<<mYBND[j+1]<<"]"<<endl;
      cout<<"pt > "<<mMinPt[itrig]<<" GeV, PFJets =  "<<pf_counter_pt[itrig][j]<<", CaloJets = "<<calo_counter_pt[itrig][j]<<endl;
      cout<<"id, PFJets = "<<pf_counter_id[itrig][j]<<", CaloJets = "<<calo_counter_id[itrig][j]<<endl;
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
InclusiveHistos::~InclusiveHistos() 
{
}

DEFINE_FWK_MODULE(InclusiveHistos);
