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
  mMaxPt     = cfg.getParameter<std::vector<double> > ("maxPt");
  mFileName  = cfg.getParameter<std::string> ("filename");
  mTreeName  = cfg.getParameter<std::string> ("treename");
  mDirName   = cfg.getParameter<std::string> ("dirname");
  mTriggers  = cfg.getParameter<std::vector<std::string> > ("triggers");
  mTrigComb  = cfg.getParameter<bool> ("trigger_and_or");    
  mUsePF     = cfg.getParameter<bool> ("usePF");
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
  mTrigIndex.clear();
  cout<<"Finding trigger mapping: "<<endl;
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
  //--------- book histos -----------------------
  TH1F *h;
  for(int i=0;i<4;i++) {
    sprintf(name,"METovSUMET_Jet%d",i+1);
    h = fs->make<TH1F>(name,name,100,0,1.0001);
    mhMETovSUMET.push_back(h);
  }
  mhJetMulti = fs->make<TH1F>("JetMultiplicity","JetMultiplicity",20,0,20);
  double aux[500];
  for(unsigned ipt=0;ipt<mPTBND.size();ipt++)
    aux[ipt] = mPTBND[ipt];
  for(unsigned iy=0;iy<mYBND.size()-1;iy++) {
    sprintf(name,"Pt_Ybin%d",iy);
    h = fs->make<TH1F>(name,name,mPTBND.size()-1,aux);
    h->Sumw2();
    mhPt.push_back(h);
    sprintf(name,"NormPt_Ybin%d",iy);
    h = fs->make<TH1F>(name,name,mPTBND.size()-1,aux);
    h->Sumw2();
    mhNormPt.push_back(h);
    sprintf(name,"TruncPt_Ybin%d",iy);
    h = fs->make<TH1F>(name,name,mPTBND.size()-1,aux);
    h->Sumw2();
    mhTruncPt.push_back(h);
    sprintf(name,"NormTruncPt_Ybin%d",iy);
    h = fs->make<TH1F>(name,name,mPTBND.size()-1,aux);
    h->Sumw2();
    mhNormTruncPt.push_back(h);
    if (mUsePF) { 
      sprintf(name,"CHF_Ybin%d",iy);
      h = fs->make<TH1F>(name,name,100,0,1.001);
      mhCHF.push_back(h);
      sprintf(name,"NHF_Ybin%d",iy);
      h = fs->make<TH1F>(name,name,100,0,1.001);
      mhNHF.push_back(h);
      sprintf(name,"PHF_Ybin%d",iy);
      h = fs->make<TH1F>(name,name,100,0,1.001);
      mhPHF.push_back(h);
    }
    else {
      sprintf(name,"EMF_Ybin%d",iy);
      h = fs->make<TH1F>(name,name,100,0,1.001);
      mhEMF.push_back(h);
      sprintf(name,"N90hits_Ybin%d",iy);
      h = fs->make<TH1F>(name,name,200,0,200);
      mhN90hits.push_back(h);
      sprintf(name,"fHPD_Ybin%d",iy);
      h = fs->make<TH1F>(name,name,100,0,1.001);
      mhfHPD.push_back(h);
      sprintf(name,"NTrkCalo_Ybin%d",iy);
      h = fs->make<TH1F>(name,name,200,0,200);
      mhNTrkCalo.push_back(h);
      sprintf(name,"fTrkVtx_Ybin%d",iy);
      h = fs->make<TH1F>(name,name,200,0,200);
      mhNTrkVtx.push_back(h);  
    }
  }
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
  int counter_hlt(0),counter_pv(0);
  int counter_y[10],counter_pt[10],counter_id[10];
  for(int iy=0;iy<10;iy++) {
    counter_y[iy] = 0;
    counter_pt[iy] = 0;
    counter_id[iy] = 0;
  }
  for(unsigned i=0;i<NEntries;i++) {
    double progress = 10.0*i/(1.0*NEntries);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;          
    mTree->GetEntry(i);
    //---------- determining the trigger decision -----------
    int prescale(1);
    bool hltPass(false);
    int tmphlt(0);
    for(unsigned itrig=0;itrig<mTrigIndex.size();itrig++) {
      int ihlt = mTrigIndex[itrig];
      if (mTrigComb) { // AND between trigger paths
        if (mEvent->fired(ihlt) > -1) {
          tmphlt *= mEvent->fired(ihlt);
          if (mTrigIndex.size() == 1)
            prescale = mEvent->preL1(ihlt) * mEvent->preHLT(ihlt);
        }
      }
      else {// OR between trigger paths
        if (mEvent->fired(ihlt) > -1) {
          tmphlt += mEvent->fired(ihlt);
          prescale = mEvent->preL1(ihlt) * mEvent->preHLT(ihlt);
        }
      }
    }
    if (tmphlt > 0)
      hltPass = true;
    //-------- cut flow ------------------------------------
    if (hltPass) {
      counter_hlt++;
      if (mEvent->evtHdr().isPVgood() == 1) {
        counter_pv++;
        if (mUsePF) {
          int nGoodJets(0);
          for(unsigned j=0;j<mEvent->nPFJets();j++) {
            int ybin = getBin(fabs(mEvent->pfjet(j).y()),mYBND);  
            if (ybin > -1) {
              counter_y[ybin]++;
              if (mEvent->pfjet(j).ptCor() >= mMinPt[ybin]) {
                counter_pt[ybin]++;
                if (mEvent->pfjet(j).id() == 1) {
                  counter_id[ybin]++;
                  nGoodJets++;
                  mhPt[ybin]->Fill((mEvent->pfjet(j)).ptCor());
                  mhNormPt[ybin]->Fill((mEvent->pfjet(j)).ptCor(),prescale);
                  mhCHF[ybin]->Fill((mEvent->pfjet(j)).chf());
                  mhNHF[ybin]->Fill((mEvent->pfjet(j)).nhf());
                  mhPHF[ybin]->Fill((mEvent->pfjet(j)).phf());
                  if (mEvent->pfjet(j).ptCor() < mMaxPt[ybin]) {
                    mhTruncPt[ybin]->Fill((mEvent->pfjet(j)).ptCor());
                    mhNormTruncPt[ybin]->Fill((mEvent->pfjet(j)).ptCor(),prescale);
                  }// cut max pt
                }// cut id
              }// cut max pt
            }// cut y  
          }// jet loop
          mhJetMulti->Fill(nGoodJets);
          for(int ii=0;ii<4;ii++) {
            if (nGoodJets > ii)
              mhMETovSUMET[ii]->Fill(mEvent->pfmet().met_o_sumet());
          }
        }// if UsePF
        else {
          int nGoodJets(0);
          for(unsigned j=0;j<mEvent->nCaloJets();j++) {
            int ybin = getBin(fabs(mEvent->calojet(j).y()),mYBND);
            if (ybin > -1) {
              counter_y[ybin]++;
              if (mEvent->calojet(j).ptCor() >= mMinPt[ybin]) {
                counter_pt[ybin]++;
                if (mEvent->calojet(j).id() == 1) {
                  counter_id[ybin]++;
                  nGoodJets++;
                  mhPt[ybin]->Fill((mEvent->calojet(j)).ptCor());
                  mhNormPt[ybin]->Fill((mEvent->calojet(j)).ptCor(),prescale);
                  mhEMF[ybin] ->Fill((mEvent->calojet(j)).emf());
                  mhN90hits[ybin]->Fill((mEvent->calojet(j)).n90hits());
                  mhfHPD[ybin]->Fill((mEvent->calojet(j)).fHPD());
                  mhNTrkCalo[ybin]->Fill((mEvent->calojet(j)).nTrkCalo());
                  mhNTrkVtx[ybin]->Fill((mEvent->calojet(j)).nTrkVtx());              
                  if (mEvent->calojet(j).ptCor() < mMaxPt[ybin]) {
                    mhTruncPt[ybin]->Fill((mEvent->calojet(j)).ptCor());
                    mhNormTruncPt[ybin]->Fill((mEvent->calojet(j)).ptCor(),prescale);
                  }// cut max pt
                }// cut id
              }// cut min pt
            }// cut y
          }// jet loop
          mhJetMulti->Fill(nGoodJets);
          for(int ii=0;ii<4;ii++) {
            if (nGoodJets > ii)
              mhMETovSUMET[ii]->Fill(mEvent->calomet().met_o_sumet());
          }
        }// if calo
      }// if pv
    }// if hlt
  }// tree loop
  cout<<"Events passing the trigger: "<<counter_hlt<<endl;
  cout<<"Events passing the PV:      "<<counter_pv<<endl;
  for(unsigned j=0;j<mYBND.size()-1;j++) {
    cout<<"--------------------------------------------"<<endl; 
    cout<<"["<<mYBND[j]<<", "<<mYBND[j+1]<<"]"<<endl;
    cout<<"Jets with pt > "<<mMinPt[j]<<" GeV: "<<counter_pt[j]<<endl;
    cout<<"Jets passing the id:   "<<counter_id[j]<<endl;
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
InclusiveHistos::~InclusiveHistos() 
{
}

DEFINE_FWK_MODULE(InclusiveHistos);
