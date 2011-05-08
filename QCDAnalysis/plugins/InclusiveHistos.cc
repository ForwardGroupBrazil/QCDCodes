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
  mMinPt     = cfg.getParameter<double> ("minPt");
  mFileName  = cfg.getParameter<std::string> ("filename");
  mTreeName  = cfg.getParameter<std::string> ("treename");
  mDirName   = cfg.getParameter<std::string> ("dirname");
  mTriggers  = cfg.getParameter<std::vector<std::string> > ("triggers");
  mUsePF     = cfg.getParameter<bool> ("usePF");
  mIsMC      = cfg.getParameter<bool> ("isMC");
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
  double auxPt[500],auxX[500];
  for(unsigned ipt=0;ipt<mPTBND.size();ipt++) {
    auxPt[ipt] = mPTBND[ipt];
    auxX[ipt]  = mPTBND[ipt]/3500.;
  }
  if (mTriggers.size() == 0) {
    for(int im=0;im<4;im++) {
      sprintf(name,"METovSUMET_%dJet",im+1);
      mhMETovSUMET[0][im] = fs->make<TH1F>(name,name,100,0,1.0001);
    }
    sprintf(name,"JetMultiplicity");
    mhJetMulti[0] = fs->make<TH1F>(name,name,20,0,20);
    for(unsigned iy=0;iy<mYBND.size()-1;iy++) {
      //-------- pt histograms ----------------------------
      sprintf(name,"Pt_Ybin%d",iy);
      mhPt[0][iy] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mhPt[0][iy]->Sumw2();
      sprintf(name,"NormPt_Ybin%d",iy);
      mhNormPt[0][iy] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mhNormPt[0][iy]->Sumw2();
      //-------- x histograms -----------------------------
      sprintf(name,"X_Ybin%d",iy);
      mhX[0][iy] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxX);
      mhX[0][iy]->Sumw2();
      sprintf(name,"NormX_Ybin%d",iy);
      mhNormX[0][iy] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxX);
      mhNormX[0][iy]->Sumw2();
      if (mUsePF) {
        sprintf(name,"CHF_Ybin%d",iy);
        mhCHF[0][iy] = fs->make<TH1F>(name,name,100,0,1.001);
        mhCHF[0][iy]->Sumw2();
        sprintf(name,"NHF_Ybin%d",iy);
        mhNHF[0][iy] = fs->make<TH1F>(name,name,100,0,1.001);
        mhNHF[0][iy]->Sumw2();
        sprintf(name,"PHF_Ybin%d",iy);
        mhPHF[0][iy] = fs->make<TH1F>(name,name,100,0,1.001);
        mhPHF[0][iy]->Sumw2();
      }
      else {
        sprintf(name,"EMF_Ybin%d",iy);
        mhEMF[0][iy] = fs->make<TH1F>(name,name,100,0,1.001);
        mhEMF[0][iy]->Sumw2();
        sprintf(name,"N90hits_Ybin%d",iy);
        mhN90hits[0][iy] = fs->make<TH1F>(name,name,200,0,200);
        mhN90hits[0][iy]->Sumw2();
        sprintf(name,"fHPD_Ybin%d",iy);
        mhfHPD[0][iy] = fs->make<TH1F>(name,name,100,0,1.001);
        mhfHPD[0][iy]->Sumw2();
        sprintf(name,"NTrkCalo_Ybin%d",iy);
        mhNTrkCalo[0][iy] = fs->make<TH1F>(name,name,200,0,200);
        mhNTrkCalo[0][iy]->Sumw2();
        sprintf(name,"fTrkVtx_Ybin%d",iy);
        mhNTrkVtx[0][iy] = fs->make<TH1F>(name,name,200,0,200);
        mhNTrkVtx[0][iy]->Sumw2();
      }
    }
  }
  for(unsigned itrig=0;itrig<mTriggers.size();itrig++) {
    for(int im=0;im<4;im++) {
      sprintf(name,"METovSUMET_%dJet_%s",im+1,mTriggers[itrig].c_str());
      mhMETovSUMET[itrig][im] = fs->make<TH1F>(name,name,100,0,1.0001);
    }
    sprintf(name,"JetMultiplicity_%s",mTriggers[itrig].c_str());
    mhJetMulti[itrig] = fs->make<TH1F>(name,name,20,0,20);
    for(unsigned iy=0;iy<mYBND.size()-1;iy++) {
      //-------- pt histograms ----------------------------
      sprintf(name,"Pt_Ybin%d_%s",iy,mTriggers[itrig].c_str());
      mhPt[itrig][iy] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mhPt[itrig][iy]->Sumw2();
      sprintf(name,"NormPt_Ybin%d_%s",iy,mTriggers[itrig].c_str());
      mhNormPt[itrig][iy] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mhNormPt[itrig][iy]->Sumw2();
      //-------- x histograms -----------------------------
      sprintf(name,"X_Ybin%d_%s",iy,mTriggers[itrig].c_str());
      mhX[itrig][iy] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxX);
      mhX[itrig][iy]->Sumw2();
      sprintf(name,"NormX_Ybin%d_%s",iy,mTriggers[itrig].c_str());
      mhNormX[itrig][iy] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxX);
      mhNormX[itrig][iy]->Sumw2();
      if (mUsePF) { 
        sprintf(name,"CHF_Ybin%d_%s",iy,mTriggers[itrig].c_str());
        mhCHF[itrig][iy] = fs->make<TH1F>(name,name,100,0,1.001);
        mhCHF[itrig][iy]->Sumw2();  
        sprintf(name,"NHF_Ybin%d_%s",iy,mTriggers[itrig].c_str());
        mhNHF[itrig][iy] = fs->make<TH1F>(name,name,100,0,1.001);
        mhNHF[itrig][iy]->Sumw2();
        sprintf(name,"PHF_Ybin%d_%s",iy,mTriggers[itrig].c_str());
        mhPHF[itrig][iy] = fs->make<TH1F>(name,name,100,0,1.001);
        mhPHF[itrig][iy]->Sumw2();
      }
      else {
        sprintf(name,"EMF_Ybin%d_%s",iy,mTriggers[itrig].c_str());
        mhEMF[itrig][iy] = fs->make<TH1F>(name,name,100,0,1.001);
        mhEMF[itrig][iy]->Sumw2();
        sprintf(name,"N90hits_Ybin%d_%s",iy,mTriggers[itrig].c_str());
        mhN90hits[itrig][iy] = fs->make<TH1F>(name,name,200,0,200);
        mhN90hits[itrig][iy]->Sumw2();
        sprintf(name,"fHPD_Ybin%d_%s",iy,mTriggers[itrig].c_str());
        mhfHPD[itrig][iy] = fs->make<TH1F>(name,name,100,0,1.001);
        mhfHPD[itrig][iy]->Sumw2();
        sprintf(name,"NTrkCalo_Ybin%d_%s",iy,mTriggers[itrig].c_str());
        mhNTrkCalo[itrig][iy] = fs->make<TH1F>(name,name,200,0,200);
        mhNTrkCalo[itrig][iy]->Sumw2();
        sprintf(name,"fTrkVtx_Ybin%d_%s",iy,mTriggers[itrig].c_str());
        mhNTrkVtx[itrig][iy] = fs->make<TH1F>(name,name,200,0,200);
        mhNTrkVtx[itrig][iy]->Sumw2();
      }
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
void InclusiveHistos::analyze(edm::Event const& evt, edm::EventSetup const& iSetup) 
{ 
  unsigned NEntries = mTree->GetEntries();
  cout<<"File: "<<mFileName<<endl;
  cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
  int decade = 0;
  int counter_hlt[30],counter_pv[30];
  int counter_y[30][10],counter_pt[30][10],counter_id[30][10];
  for(unsigned itrig=0;itrig<mTrigIndex.size();itrig++) {
    counter_hlt[itrig] = 0;
    counter_pv[itrig]  = 0;
    for(int iy=0;iy<10;iy++) {
      counter_y[itrig][iy] = 0;
      counter_pt[itrig][iy] = 0;
      counter_id[itrig][iy] = 0;
    }
  }
  if (mTriggers.size() == 0) {
    counter_pv[0]  = 0;
    for(int iy=0;iy<10;iy++) {
      counter_y[0][iy] = 0;
      counter_pt[0][iy] = 0;
      counter_id[0][iy] = 0;
    }
  }
  //---------- loop over the events ----------------------
  for(unsigned i=0;i<NEntries;i++) {
    double progress = 10.0*i/(1.0*NEntries);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;          
    mTree->GetEntry(i);
    double wt(1.0);
    if (mIsMC)
      wt = mEvent->evtHdr().weight(); 
    if (mTriggers.size() == 0) {
      int prescale(1.0);
      if (mEvent->evtHdr().isPVgood() == 1) {
        counter_pv[0]++;
        if (mUsePF) {
          int nGoodJets(0);
          for(unsigned j=0;j<mEvent->nPFJets();j++) {
            int ybin = getBin(fabs(mEvent->pfjet(j).y()),mYBND);
            if (ybin > -1) {
              counter_y[0][ybin]++;
              if (mEvent->pfjet(j).ptCor() >= mMinPt) {
                counter_pt[0][ybin]++;
                if (mEvent->pfjet(j).id() == 1) {
                  counter_id[0][ybin]++;
                  nGoodJets++;
                  mhPt[0][ybin]->Fill((mEvent->pfjet(j)).ptCor(),wt);
                  mhNormPt[0][ybin]->Fill((mEvent->pfjet(j)).ptCor(),prescale);
                  double x = mEvent->pfjet(j).ptCor()*cosh(mEvent->pfjet(j).y())/3500;
                  mhX[0][ybin]->Fill(x,wt);
                  mhNormX[0][ybin]->Fill(x,prescale);
                  mhCHF[0][ybin]->Fill((mEvent->pfjet(j)).chf(),wt);
                  mhNHF[0][ybin]->Fill((mEvent->pfjet(j)).nhf(),wt);
                  mhPHF[0][ybin]->Fill((mEvent->pfjet(j)).phf(),wt);
                }// cut id
              }// cut max pt
            }// cut y  
          }// jet loop
          mhJetMulti[0]->Fill(nGoodJets,wt);
          for(int ii=0;ii<4;ii++) {
            if (nGoodJets > ii)
              mhMETovSUMET[0][ii]->Fill(mEvent->pfmet().met_o_sumet(),wt);
          }
        }// if UsePF
        else {
          int nGoodJets(0);
          for(unsigned j=0;j<mEvent->nCaloJets();j++) {
            int ybin = getBin(fabs(mEvent->calojet(j).y()),mYBND);
            if (ybin > -1) {
              counter_y[0][ybin]++;
              if (mEvent->calojet(j).ptCor() >= mMinPt) {
                counter_pt[0][ybin]++;
                if (mEvent->calojet(j).id() == 1) {
                  counter_id[0][ybin]++;
                  nGoodJets++;
                  mhPt[0][ybin]->Fill((mEvent->calojet(j)).ptCor(),wt);
                  mhNormPt[0][ybin]->Fill((mEvent->calojet(j)).ptCor(),prescale);
                  double x = mEvent->calojet(j).ptCor()*cosh(mEvent->calojet(j).y())/3500;
                  mhX[0][ybin]->Fill(x,wt);
                  mhNormX[0][ybin]->Fill(x,prescale);
                  mhEMF[0][ybin] ->Fill((mEvent->calojet(j)).emf(),wt);
                  mhN90hits[0][ybin]->Fill((mEvent->calojet(j)).n90hits(),wt);
                  mhfHPD[0][ybin]->Fill((mEvent->calojet(j)).fHPD(),wt);
                  mhNTrkCalo[0][ybin]->Fill((mEvent->calojet(j)).nTrkCalo(),wt);
                  mhNTrkVtx[0][ybin]->Fill((mEvent->calojet(j)).nTrkVtx(),wt);
                }// cut id
              }// cut min pt
            }// cut y
          }// jet loop
          mhJetMulti[0]->Fill(nGoodJets,wt);
          for(int ii=0;ii<4;ii++) {
            if (nGoodJets > ii)
              mhMETovSUMET[0][ii]->Fill(mEvent->calomet().met_o_sumet(),wt);
          }
        }// if calo
      }// if pv
    }
    //---------- loop over the triggers ----------------------
    for(unsigned itrig=0;itrig<mTrigIndex.size();itrig++) {
      int prescale(1);
      bool hltPass(false);
      int ihlt = mTrigIndex[itrig];
      if (mEvent->fired(ihlt) > 0) {
        hltPass = true;
        prescale = mEvent->preL1(ihlt) * mEvent->preHLT(ihlt);
      }
      //-------- cut flow ------------------------------------
      if (hltPass) {
        counter_hlt[itrig]++;
        if (mEvent->evtHdr().isPVgood() == 1) {
          counter_pv[itrig]++;
          if (mUsePF) {
            int nGoodJets(0);
            for(unsigned j=0;j<mEvent->nPFJets();j++) {
              int ybin = getBin(fabs(mEvent->pfjet(j).y()),mYBND);  
              if (ybin > -1) {
                counter_y[itrig][ybin]++;
                if (mEvent->pfjet(j).ptCor() >= mMinPt) {
                  counter_pt[itrig][ybin]++;
                  if (mEvent->pfjet(j).id() == 1) {
                    counter_id[itrig][ybin]++;
                    nGoodJets++;
                    mhPt[itrig][ybin]->Fill((mEvent->pfjet(j)).ptCor(),wt);
                    mhNormPt[itrig][ybin]->Fill((mEvent->pfjet(j)).ptCor(),prescale);
                    double x = mEvent->pfjet(j).ptCor()*cosh(mEvent->pfjet(j).y())/3500;
                    mhX[itrig][ybin]->Fill(x,wt);
                    mhNormX[itrig][ybin]->Fill(x,prescale);
                    mhCHF[itrig][ybin]->Fill((mEvent->pfjet(j)).chf(),wt);
                    mhNHF[itrig][ybin]->Fill((mEvent->pfjet(j)).nhf(),wt);
                    mhPHF[itrig][ybin]->Fill((mEvent->pfjet(j)).phf(),wt);
                  }// cut id
                }// cut max pt
              }// cut y  
            }// jet loop
            mhJetMulti[itrig]->Fill(nGoodJets,wt);
            for(int ii=0;ii<4;ii++) {
              if (nGoodJets > ii)
                mhMETovSUMET[itrig][ii]->Fill(mEvent->pfmet().met_o_sumet(),wt);
            }
          }// if UsePF
          else {
            int nGoodJets(0);
            for(unsigned j=0;j<mEvent->nCaloJets();j++) {
              int ybin = getBin(fabs(mEvent->calojet(j).y()),mYBND);
              if (ybin > -1) {
                counter_y[itrig][ybin]++;
                if (mEvent->calojet(j).ptCor() >= mMinPt) {
                  counter_pt[itrig][ybin]++;
                  if (mEvent->calojet(j).id() == 1) {
                    counter_id[itrig][ybin]++;
                    nGoodJets++;
                    mhPt[itrig][ybin]->Fill((mEvent->calojet(j)).ptCor(),wt);
                    mhNormPt[itrig][ybin]->Fill((mEvent->calojet(j)).ptCor(),prescale);
                    double x = mEvent->calojet(j).ptCor()*cosh(mEvent->calojet(j).y())/3500;
                    mhX[itrig][ybin]->Fill(x,wt);
                    mhNormX[itrig][ybin]->Fill(x,prescale);
                    mhEMF[itrig][ybin] ->Fill((mEvent->calojet(j)).emf(),wt);
                    mhN90hits[itrig][ybin]->Fill((mEvent->calojet(j)).n90hits(),wt);
                    mhfHPD[itrig][ybin]->Fill((mEvent->calojet(j)).fHPD(),wt);
                    mhNTrkCalo[itrig][ybin]->Fill((mEvent->calojet(j)).nTrkCalo(),wt);
                    mhNTrkVtx[itrig][ybin]->Fill((mEvent->calojet(j)).nTrkVtx(),wt);              
                  }// cut id
                }// cut min pt
              }// cut y
            }// jet loop
            mhJetMulti[itrig]->Fill(nGoodJets,wt);
            for(int ii=0;ii<4;ii++) {
              if (nGoodJets > ii)
                mhMETovSUMET[itrig][ii]->Fill(mEvent->calomet().met_o_sumet(),wt);
            }
          }// if calo
        }// if pv
      }// if hlt
    }// trigger loop
  }// tree loop
  for(unsigned itrig=0;itrig<mTrigIndex.size();itrig++) {
    cout<<"*********************************************"<<endl;
    cout<<"Trigger path: "              <<mTriggers[itrig]<<endl;
    cout<<"Events passing the trigger: "<<counter_hlt[itrig]<<endl;
    cout<<"Events passing the PV:      "<<counter_pv[itrig]<<endl;
    for(unsigned j=0;j<mYBND.size()-1;j++) {
      cout<<"--------------------------------------------"<<endl; 
      cout<<"["<<mYBND[j]<<", "<<mYBND[j+1]<<"]"<<endl;
      cout<<"Jets with pt > "<<mMinPt<<" GeV: "<<counter_pt[itrig][j]<<endl;
      cout<<"Jets passing the id:   "<<counter_id[itrig][j]<<endl;
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
InclusiveHistos::~InclusiveHistos() 
{
}

DEFINE_FWK_MODULE(InclusiveHistos);
