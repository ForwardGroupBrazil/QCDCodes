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

#include "KKousour/QCDAnalysis/plugins/DijetMassHistos.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

DijetMassHistos::DijetMassHistos(edm::ParameterSet const& cfg) 
{
  mYBND      = cfg.getParameter<std::vector<double> > ("yBnd");
  mMASSBND   = cfg.getParameter<std::vector<double> > ("massBnd");
  mMinMass   = cfg.getParameter<std::vector<double> > ("minMass");
  mMaxMass   = cfg.getParameter<std::vector<double> > ("maxMass");
  mMinPt1    = cfg.getParameter<double> ("minPt1");
  mMinPt2    = cfg.getParameter<double> ("minPt2");
  mFileName  = cfg.getParameter<std::string> ("filename");
  mTreeName  = cfg.getParameter<std::string> ("treename");
  mUsePF     = cfg.getParameter<bool> ("usePF");
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetMassHistos::beginJob() 
{
  mInf = TFile::Open(mFileName.c_str());
  mTree = (TTree*)mInf->Get(mTreeName.c_str());
  mEvent = new QCDEvent();
  TBranch *branch = mTree->GetBranch("event");
  branch->SetAddress(&mEvent);
  //--------- book histos -----------------------
  char name[1000];
  TH1F *h;
  double aux[500];
  for(unsigned imass=0;imass<mMASSBND.size();imass++)
    aux[imass] = mMASSBND[imass];
  for(unsigned iy=0;iy<mYBND.size()-1;iy++) {
    sprintf(name,"METovSUMET_Ybin%d",iy);
    h = fs->make<TH1F>("METovSUMET","METovSUMET",100,0,1.0001); 
    mhMETovSUMET.push_back(h);
    sprintf(name,"Mass_Ybin%d",iy);
    h = fs->make<TH1F>(name,name,mMASSBND.size()-1,aux);
    h->Sumw2();
    mhM.push_back(h);
    sprintf(name,"NormMass_Ybin%d",iy);
    h = fs->make<TH1F>(name,name,mMASSBND.size()-1,aux);
    h->Sumw2();
    mhNormM.push_back(h);
    sprintf(name,"TruncMass_Ybin%d",iy);
    h = fs->make<TH1F>(name,name,mMASSBND.size()-1,aux);
    h->Sumw2();
    mhTruncM.push_back(h);
    sprintf(name,"NormTruncMass_Ybin%d",iy);
    h = fs->make<TH1F>(name,name,mMASSBND.size()-1,aux);
    h->Sumw2();
    mhNormTruncM.push_back(h);
    sprintf(name,"Pt_Ybin%d",iy);
    h = fs->make<TH1F>(name,name,3500,0,3500);
    h->Sumw2();
    mhPt.push_back(h);
    sprintf(name,"Y_Ybin%d",iy);
    h = fs->make<TH1F>(name,name,200,-5,5);
    h->Sumw2();
    mhY.push_back(h);
    sprintf(name,"Ymax_Ybin%d",iy);
    h = fs->make<TH1F>(name,name,200,0,5);
    h->Sumw2();
    mhYmax.push_back(h);
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
void DijetMassHistos::endJob() 
{
  mInf->Close();
}
//////////////////////////////////////////////////////////////////////////////////////////
int DijetMassHistos::getBin(double x, const std::vector<double>& boundaries)
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
void DijetMassHistos::analyze(edm::Event const& evt, edm::EventSetup const& iSetup) 
{ 
  unsigned NEntries = mTree->GetEntries();
  cout<<"File: "<<mFileName<<endl;
  cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
  int decade = 0;
  for(unsigned i=0;i<NEntries;i++) {
    double progress = 10.0*i/(1.0*NEntries);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;          
    mTree->GetEntry(i);
    if (mUsePF) {
      if (mEvent->evtHdr().isPVgood() == 1 && mEvent->nPFJets() > 1 ) {   
        int prescale = mEvent->evtHdr().preL1() * mEvent->evtHdr().preHLT();
        double ymax = TMath::Max(fabs(mEvent->pfjet(0).y()),fabs(mEvent->pfjet(1).y()));
        int ybin = getBin(ymax,mYBND);
        bool cut1 = (mEvent->pfjet(0).id() == 1 && mEvent->pfjet(0).ptCor() > mMinPt1);
        bool cut2 = (mEvent->pfjet(1).id() == 1 && mEvent->pfjet(1).ptCor() > mMinPt2);
        if (cut1 && cut2 && ybin > -1) {
          double mjj = mEvent->pfmjjcor(0);
          if (mjj >= mMinMass[ybin]) {
            mhM[ybin]->Fill(mjj);
            mhNormM[ybin]->Fill(mjj,prescale);
            if (mjj < mMaxMass[ybin]) {
              mhMETovSUMET[ybin]->Fill(mEvent->pfmet().met_o_sumet());
              mhTruncM[ybin]->Fill(mjj);
              mhNormTruncM[ybin]->Fill(mjj,prescale);
              mhYmax[ybin]->Fill(ymax);             
              for(unsigned j=0;j<2;j++) {
                mhPt[ybin]->Fill((mEvent->pfjet(j)).ptCor());
                mhY[ybin]->Fill((mEvent->pfjet(j)).y());
                mhCHF[ybin]->Fill((mEvent->pfjet(j)).chf());
                mhNHF[ybin]->Fill((mEvent->pfjet(j)).nhf());
                mhPHF[ybin]->Fill((mEvent->pfjet(j)).phf());
              }
            }
          }
        }
      }  
    }
    else {
      if (mEvent->evtHdr().isPVgood() == 1 && mEvent->nCaloJets() > 1 ) {
        int prescale = mEvent->evtHdr().preL1() * mEvent->evtHdr().preHLT();
        double ymax = TMath::Max(fabs(mEvent->calojet(0).y()),fabs(mEvent->calojet(1).y()));
        int ybin = getBin(ymax,mYBND);
        bool cut1 = (mEvent->calojet(0).id() == 1 && mEvent->calojet(0).ptCor() > mMinPt1);
        bool cut2 = (mEvent->calojet(1).id() == 1 && mEvent->calojet(1).ptCor() > mMinPt2);
        if (cut1 && cut2 && ybin > -1) {
          double mjj = mEvent->calomjjcor(0);
          if (mjj >= mMinMass[ybin]) {
            mhM[ybin]->Fill(mjj);
            mhNormM[ybin]->Fill(mjj,prescale);
            if (mjj < mMaxMass[ybin]) {
              mhTruncM[ybin]->Fill(mjj);
              mhNormTruncM[ybin]->Fill(mjj,prescale);
              mhYmax[ybin]->Fill(ymax);
              mhMETovSUMET[ybin]->Fill(mEvent->calomet().met_o_sumet());
              for(unsigned j=0;j<2;j++) {
                mhPt[ybin]->Fill((mEvent->calojet(j)).ptCor());
                mhY[ybin]->Fill((mEvent->calojet(j)).y());
                mhEMF[ybin] ->Fill((mEvent->calojet(j)).emf());
                mhN90hits[ybin]->Fill((mEvent->calojet(j)).n90hits());
                mhfHPD[ybin]->Fill((mEvent->calojet(j)).fHPD());
                mhNTrkCalo[ybin]->Fill((mEvent->calojet(j)).nTrkCalo());
                mhNTrkVtx[ybin]->Fill((mEvent->calojet(j)).nTrkVtx());
              }
            }
          }
        }
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////

DEFINE_FWK_MODULE(DijetMassHistos);
