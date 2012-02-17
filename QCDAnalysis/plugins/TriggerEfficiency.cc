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

#include "KKousour/QCDAnalysis/plugins/TriggerEfficiency.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

TriggerEfficiency::TriggerEfficiency(edm::ParameterSet const& cfg) 
{
  mYBND       = cfg.getParameter<std::vector<double> >      ("yBnd");
  mPTBND      = cfg.getParameter<std::vector<double> >      ("ptBnd");
  mMinPt      = cfg.getParameter<double>                    ("minPt");
  mL1Pt       = cfg.getParameter<std::vector<double> >      ("L1Pt");
  mHLTPt      = cfg.getParameter<std::vector<double> >      ("HLTPt");
  mFileNames  = cfg.getParameter<std::vector<std::string> > ("filenames");
  mTreeName   = cfg.getParameter<std::string>               ("treename");
  mDirName    = cfg.getParameter<std::string>               ("dirname");
  mRefTrigger = cfg.getParameter<std::vector<std::string> > ("refTrigger");
  mNEvents    = cfg.getParameter<int>                       ("nEvents"); 
  mJetID      = cfg.getParameter<int>                       ("jetID");
  mHCALNoise  = cfg.getParameter<int>                       ("hcalNoiseFilter");
  if (mL1Pt.size() != mHLTPt.size() || mL1Pt.size() != mRefTrigger.size())
    throw cms::Exception("TriggerEfficiency: ")<<"inconsistent trigger settings";
  cout<<"Number of triggers to be emulated: "<<mL1Pt.size()<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void TriggerEfficiency::beginJob() 
{
  mInf = TFile::Open(mFileNames[0].c_str());
  mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());
  mEvent = new QCDEvent();
  char name[1000];
  //--------- trigger mapping -------------------
  TH1F *hTrigNames = (TH1F*)mDir->Get("TriggerNames");
  cout<<"Finding trigger mapping: "<<endl;
  for(unsigned itrig=0;itrig<mRefTrigger.size();itrig++) {
    int index(-1);
    for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
      string ss = hTrigNames->GetXaxis()->GetBinLabel(ibin+1);
      if (ss == mRefTrigger[itrig]) {
        index = ibin;
        continue;
      }
    }
    if (index > -1)
      mRefTrigIndex.push_back(index); 
    else {
      throw cms::Exception("TriggerEfficiency: ")<<"The requested trigger ("<<mRefTrigger[itrig]<<") is not found ";
    }
    cout<<mRefTrigger[itrig]<<" --> "<<mRefTrigIndex[itrig]<<endl; 
  }
  //--------- book histos ---------------------------------
  double vpt[100];
  for(unsigned ipt=0;ipt<mPTBND.size();ipt++)
    vpt[ipt] = mPTBND[ipt];
  TFileDirectory mPFDir   = fs->mkdir(mDirName+"pf");
  TFileDirectory mCaloDir = fs->mkdir(mDirName+"calo");
  for(unsigned itrig=0;itrig<mRefTrigger.size();itrig++) {
    for(unsigned iy=0;iy<mYBND.size()-1;iy++) {
      //-------- pt histograms ----------------------------
      sprintf(name,"Pt_Ybin%d_%s_hlt%d",iy,mRefTrigger[itrig].c_str(),(int)mHLTPt[itrig]);
      mPFPt[itrig][iy] = mPFDir.make<TH1F>(name,name,mPTBND.size()-1,vpt);
      mPFPt[itrig][iy]->Sumw2();
      sprintf(name,"RefPt_Ybin%d_%s_hlt%d",iy,mRefTrigger[itrig].c_str(),(int)mHLTPt[itrig]);
      mPFRefPt[itrig][iy] = mPFDir.make<TH1F>(name,name,mPTBND.size()-1,vpt);
      mPFRefPt[itrig][iy]->Sumw2();
      sprintf(name,"Pt_Ybin%d_%s_hlt%d",iy,mRefTrigger[itrig].c_str(),(int)mHLTPt[itrig]);
      mCaloPt[itrig][iy] = mCaloDir.make<TH1F>(name,name,mPTBND.size()-1,vpt);
      mCaloPt[itrig][iy]->Sumw2();
      sprintf(name,"RefPt_Ybin%d_%s_hlt%d",iy,mRefTrigger[itrig].c_str(),(int)mHLTPt[itrig]);
      mCaloRefPt[itrig][iy] = mCaloDir.make<TH1F>(name,name,mPTBND.size()-1,vpt);
      mCaloRefPt[itrig][iy]->Sumw2();
    }// y loop
  }// trigger loop
}
//////////////////////////////////////////////////////////////////////////////////////////
void TriggerEfficiency::endJob() 
{
  mInf->Close();
}
//////////////////////////////////////////////////////////////////////////////////////////
int TriggerEfficiency::getBin(double x, const std::vector<double>& boundaries)
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
void TriggerEfficiency::analyze(edm::Event const& evt, edm::EventSetup const& iSetup) 
{
  for(unsigned iFile=0;iFile<mFileNames.size();iFile++) {
    TFile *inf = TFile::Open(mFileNames[iFile].c_str());
    TTree *mTree = (TTree*)inf->Get((mDirName+"/"+mTreeName).c_str()); 
    TBranch *branch = mTree->GetBranch("events");
    branch->SetAddress(&mEvent);
    unsigned NEntries = mTree->GetEntries();
    cout<<"File: "<<mFileNames[iFile]<<endl;
    cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
    int decade = 0;
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
      //------ loop over the triggers ----------------------
      for(unsigned itrig=0;itrig<mRefTrigger.size();itrig++) {
        int ihlt = mRefTrigIndex[itrig]; 
        bool hltPass(false);
        if (mEvent->fired(ihlt) > 0) {
          hltPass = true;
        }
        //-------- cut flow ------------------------------------
        bool emulateL1(false),emulateHLT(false);
        if (hltPass) {
          if (mEvent->evtHdr().isPVgood() == 1) {
            bool cut_hcalNoise(true);
            cut_hcalNoise = mEvent->evtHdr().hcalNoise();
            if (cut_hcalNoise) {
              //-------- loop over the L1 objects ------------------
              for(unsigned j=0;j<mEvent->nL1Obj(ihlt);j++) {
                if (mEvent->l1obj(ihlt,j).pt() >= mL1Pt[itrig]) {
                  emulateL1 = true;
                  continue;
                }
              }
              //-------- loop over the HLT objects ------------------
              for(unsigned j=0;j<mEvent->nHLTObj(ihlt);j++) {
                if (mEvent->hltobj(ihlt,j).pt() >= mHLTPt[itrig]) {
                  emulateHLT = true;
                continue;
                }
              }
              //-------- loop over the PFJets ------------------
              for(unsigned j=0;j<mEvent->nPFJets();j++) {
                int ybin = getBin(fabs(mEvent->pfjet(j).y()),mYBND);  
                if (ybin > -1) {
                  if (mEvent->pfjet(j).ptCor() >= mMinPt) {
                    bool cutID(true);
                    if (mJetID == 1)
                      cutID = mEvent->pfjet(j).looseID();
                    if (mJetID == 2)
                      cutID = mEvent->pfjet(j).tightID();  
                    if (cutID) {
                      mPFRefPt[itrig][ybin]->Fill((mEvent->pfjet(j)).ptCor());
                      if (emulateL1 && emulateHLT) {
                        mPFPt[itrig][ybin]->Fill((mEvent->pfjet(j)).ptCor());
                      }
                    }// cut id
                  }// cut max pt
                }// cut y  
              }// pfjet loop
              //-------- loop over the CaloJets ------------------
              for(unsigned j=0;j<mEvent->nCaloJets();j++) {
                int ybin = getBin(fabs(mEvent->calojet(j).y()),mYBND);
                if (ybin > -1) {
                  if (mEvent->calojet(j).ptCor() >= mMinPt) {
                    bool cutID(true);
                    if (mJetID == 1)
                      cutID = mEvent->calojet(j).looseID();
                    if (mJetID == 2)
                      cutID = mEvent->calojet(j).tightID();
                    if (cutID) {
                      mCaloRefPt[itrig][ybin]->Fill((mEvent->calojet(j)).ptCor());
                      if (emulateL1 && emulateHLT) {
                        mCaloPt[itrig][ybin]->Fill((mEvent->calojet(j)).ptCor());
                      }
                    }// cut id
                  }// cut min pt
                }// cut y
              }// calojet loop
            }// if hcal noise
          }// if pv
        }// if hlt
      }// trigger loop
    }// tree loop
  }// file loop
}
//////////////////////////////////////////////////////////////////////////////////////////
TriggerEfficiency::~TriggerEfficiency() 
{
}

DEFINE_FWK_MODULE(TriggerEfficiency);
