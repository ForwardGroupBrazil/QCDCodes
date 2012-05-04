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

#include "KKousour/QCDAnalysis/plugins/InclusiveJESHistos.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

InclusiveJESHistos::InclusiveJESHistos(edm::ParameterSet const& cfg) 
{
  mYBND           = cfg.getParameter<std::vector<double> >      ("yBnd");
  mPTBND          = cfg.getParameter<std::vector<double> >      ("ptBnd");
  mMinPt          = cfg.getParameter<std::vector<double> >      ("minPt");
  mMaxMETovSumET  = cfg.getParameter<double>                    ("maxMETovSumET"); 
  mFileNames      = cfg.getParameter<std::vector<std::string> > ("filenames");
  mTreeName       = cfg.getParameter<std::string>               ("treename");
  mDirName        = cfg.getParameter<std::string>               ("dirname");
  mTriggers       = cfg.getParameter<std::vector<std::string> > ("triggers");
  mNEvents        = cfg.getParameter<int>                       ("nEvents");
  mJECUncSrcNames = cfg.getParameter<std::vector<std::string> > ("jecUncSrcNames"); 
  if (mTriggers.size() > 0) {
    if (mMinPt.size() != mTriggers.size()+1) {
      throw cms::Exception("InclusiveJESHistos: ")<<" Number of pt thresholds must be equal to the number of triggers\n";
    }
  }
  else {
    if (mMinPt.size() != 1) {
      throw cms::Exception("InclusiveJESHistos: ")<<" Number of pt thresholds must be equal to 1 when no triggers are defined\n";
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void InclusiveJESHistos::beginJob() 
{
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
  double auxPt[500];
  int Ntrig = (int)mTriggers.size();
  for(unsigned ipt=0;ipt<mPTBND.size();ipt++) {
    auxPt[ipt] = mPTBND[ipt];
  }
  for(int itrig=0;itrig<TMath::Max(Ntrig,1);itrig++) {
    string ss("");
    if (Ntrig > 0) {
      ss = "_"+mTriggers[itrig].erase(mTriggers[itrig].find("v")-1,mTriggers[itrig].find("v"));
    }
    for(unsigned iy=0;iy<mYBND.size()-1;iy++) { 
      sprintf(name,"Pt_Ybin%d%s",iy,ss.c_str());
      mPt[itrig][iy] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mPt[itrig][iy]->Sumw2();
      for(unsigned isrc=0;isrc<mJECUncSrcNames.size();isrc++) {  
        sprintf(name,"PtUp_Ybin%d%s_%s",iy,ss.c_str(),mJECUncSrcNames[isrc].c_str());
        mPtUp[itrig][iy][isrc] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxPt);
        mPtUp[itrig][iy][isrc]->Sumw2();
        sprintf(name,"PtDown_Ybin%d%s_%s",iy,ss.c_str(),mJECUncSrcNames[isrc].c_str());
        mPtDown[itrig][iy][isrc] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxPt);
        mPtDown[itrig][iy][isrc]->Sumw2();
      }
    }// y loop
  }// trigger loop
  for(unsigned iy=0;iy<mYBND.size()-1;iy++) { 
    sprintf(name,"Pt_Ybin%d",iy);
    mPtTotal[iy] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxPt);
    mPtTotal[iy]->Sumw2();
    for(unsigned isrc=0;isrc<mJECUncSrcNames.size();isrc++) {  
      sprintf(name,"PtUp_Ybin%d_%s",iy,mJECUncSrcNames[isrc].c_str());
      mPtUpTotal[iy][isrc] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mPtUpTotal[iy][isrc]->Sumw2();
      sprintf(name,"PtDown_Ybin%d_%s",iy,mJECUncSrcNames[isrc].c_str());
      mPtDownTotal[iy][isrc] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mPtDownTotal[iy][isrc]->Sumw2();
      sprintf(name,"PtRatioUp_Ybin%d_%s",iy,mJECUncSrcNames[isrc].c_str());
      mPtRatioUpTotal[iy][isrc] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mPtRatioUpTotal[iy][isrc]->Sumw2();
      sprintf(name,"PtRatioDown_Ybin%d_%s",iy,mJECUncSrcNames[isrc].c_str());
      mPtRatioDownTotal[iy][isrc] = fs->make<TH1F>(name,name,mPTBND.size()-1,auxPt);
      mPtRatioDownTotal[iy][isrc]->Sumw2();
    }
  }// y loop
}
//////////////////////////////////////////////////////////////////////////////////////////
void InclusiveJESHistos::endJob() 
{
  int Ntrig = (int)mTriggers.size();
  for(unsigned iy=0;iy<mYBND.size()-1;iy++) {
    for(int itrig=0;itrig<TMath::Max(Ntrig,1);itrig++) {
      mPtTotal[iy]->Add(mPt[itrig][iy]); 
      for(unsigned isrc=0;isrc<mJECUncSrcNames.size();isrc++) { 
        mPtUpTotal[iy][isrc]->Add(mPtUp[itrig][iy][isrc]);
        mPtDownTotal[iy][isrc]->Add(mPtDown[itrig][iy][isrc]);
      }
    }
    for(unsigned isrc=0;isrc<mJECUncSrcNames.size();isrc++) { 
      mPtRatioUpTotal[iy][isrc]->Add(mPtUpTotal[iy][isrc]);
      mPtRatioUpTotal[iy][isrc]->Divide(mPtTotal[iy]);
      mPtRatioDownTotal[iy][isrc]->Add(mPtDownTotal[iy][isrc]);
      mPtRatioDownTotal[iy][isrc]->Divide(mPtTotal[iy]);
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
int InclusiveJESHistos::getBin(double x, const std::vector<double>& boundaries)
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
void InclusiveJESHistos::analyze(edm::Event const& evt, edm::EventSetup const& iSetup) 
{ 
  int decade = 0;
  int Ntrig = (int)mTriggers.size();
  //---------- loop over the files -----------------------
  for(unsigned iFile=0;iFile<mFileNames.size();iFile++) {
    TFile *inf = TFile::Open(mFileNames[iFile].c_str());  
    TTree *mTree = (TTree*)inf->Get((mDirName+"/"+mTreeName).c_str());
    unsigned long NEntries = mTree->GetEntries();
    cout<<"Adding file: "<<mFileNames[iFile]<<endl;
    cout<<"Reading TREE: "<<NEntries<<" events"<<endl;   
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
        //-------- cut flow ------------------------------------
        if (hltPass) {
          if (mEvent->pfmet().met_o_sumet() < mMaxMETovSumET) {
            if (mEvent->evtHdr().pfRho() < 100) {
              for(unsigned j=0;j<mEvent->nPFJets();j++) {
                double y = mEvent->pfjet(j).y();
                int ybin = getBin(fabs(y),mYBND);  
                if (ybin < 0) continue;
                bool cutID=(mEvent->pfjet(j).tightID() && (mEvent->pfjet(j).elf() < 0.9) && (mEvent->pfjet(j).muf() < 0.9) && (mEvent->pfjet(j).nhf() < 0.9));
                if (!cutID) continue;
                double ptCor = mEvent->pfjet(j).ptCor();
                if (ptCor >= mMinPt[itrig] && ptCor < mMinPt[itrig+1]) {
                  mPt[itrig][ybin]->Fill(ptCor,wt); 
                }
                for(unsigned isrc=0;isrc<mJECUncSrcNames.size();isrc++) {
                  double unc = mEvent->pfjet(j).uncSrc(isrc);
                  double ptUp = (1+unc)*ptCor;
                  double ptDown = (1-unc)*ptCor;   
                  if (ptUp >= mMinPt[itrig] && ptUp < mMinPt[itrig+1]) {
                    mPtUp[itrig][ybin][isrc]->Fill(ptUp,wt);
                  }
                  if (ptDown >= mMinPt[itrig] && ptDown < mMinPt[itrig+1]) {
                    mPtDown[itrig][ybin][isrc]->Fill(ptDown,wt);
                  }
                }
              }// pfjet loop  
            }// pf rho cut
          }// pf met cut  
        }// if hlt
      }// trigger loop
    }// tree loop
  }// file loop
}
//////////////////////////////////////////////////////////////////////////////////////////
InclusiveJESHistos::~InclusiveJESHistos() 
{
}

DEFINE_FWK_MODULE(InclusiveJESHistos);
