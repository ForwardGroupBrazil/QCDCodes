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
  mMinPFPt   = cfg.getParameter<double> ("minPFPt");
  mMinCaloPt = cfg.getParameter<double> ("minCaloPt");
  mFileName  = cfg.getParameter<std::string> ("filename");
  mTreeName  = cfg.getParameter<std::string> ("treename");
}
//////////////////////////////////////////////////////////////////////////////////////////
void InclusiveHistos::beginJob() 
{
  mInf = TFile::Open(mFileName.c_str());
  mTree = (TTree*)mInf->Get(mTreeName.c_str());
  mEvent = new QCDEvent();
  TBranch *branch = mTree->GetBranch("event");
  branch->SetAddress(&mEvent);
  //--------- book histos -----------------------
  char name[1000];
  mhMETovSUMET = fs->make<TH1F>("METovSUMET","METovSUMET",100,0,1.0001);
  TH1F *h;
  double aux[500];
  for(unsigned ipt=0;ipt<mPTBND.size();ipt++)
    aux[ipt] = mPTBND[ipt];
  for(unsigned iy=0;iy<mYBND.size()-1;iy++) {
    sprintf(name,"CorPt_Ybin_%d",iy);
    h = fs->make<TH1F>(name,name,mPTBND.size()-1,aux);
    h->Sumw2();
    mhPt.push_back(h);
    sprintf(name,"CHF_Ybin_%d",iy);
    h = fs->make<TH1F>(name,name,100,0,1.001);
    mhCHF.push_back(h);
    sprintf(name,"NHF_Ybin_%d",iy);
    h = fs->make<TH1F>(name,name,100,0,1.001);
    mhNHF.push_back(h);
    sprintf(name,"PHF_Ybin_%d",iy);
    h = fs->make<TH1F>(name,name,100,0,1.001);
    mhPHF.push_back(h);
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
  cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
  int decade = 0;
  for(unsigned i=0;i<NEntries;i++) {
    double progress = 10.0*i/(1.0*NEntries);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;          
    mTree->GetEntry(i);
    if (mEvent->evtHdr().isPVgood() == 1) {
      mhMETovSUMET->Fill(mEvent->pfmet().met_o_sumet());
      for(unsigned j=0;j<mEvent->nPFJets();j++) {
        int ybin = getBin(fabs(mEvent->pfjet(j).y()),mYBND);  
        if (mEvent->pfjet(j).id() == 1 && mEvent->pfjet(j).ptCor() >= mMinPFPt && ybin > -1) {
          mhPt[ybin] ->Fill((mEvent->pfjet(j)).ptCor());
          mhCHF[ybin]->Fill((mEvent->pfjet(j)).chf());
          mhNHF[ybin]->Fill((mEvent->pfjet(j)).nhf());
          mhPHF[ybin]->Fill((mEvent->pfjet(j)).phf());
        }  
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
InclusiveHistos::~InclusiveHistos() 
{
}

DEFINE_FWK_MODULE(InclusiveHistos);
