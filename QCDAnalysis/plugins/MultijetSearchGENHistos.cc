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

#include "KKousour/QCDAnalysis/plugins/MultijetSearchGENHistos.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

MultijetSearchGENHistos::MultijetSearchGENHistos(edm::ParameterSet const& cfg) 
{
  mFileName  = cfg.getParameter<std::string>               ("filename");
  mTreeName  = cfg.getParameter<std::string>               ("treename");
  mDirName   = cfg.getParameter<std::string>               ("dirname");
  mRank      = cfg.getParameter<int>                       ("rank");
  mNEvents   = cfg.getParameter<int>                       ("nEvents");
  mScale4J   = cfg.getParameter<double>                    ("scale4j");
  mOffset4J  = cfg.getParameter<double>                    ("offset4j"); 
  mMaxEta    = cfg.getParameter<double>                    ("maxEta");
  mMinHT     = cfg.getParameter<double>                    ("minHT");
  mMinPt     = cfg.getParameter<std::vector<double> >      ("minPt");
  vector<double> vLumi(0),vBnd(0);
  mPtHatLumi = cfg.getUntrackedParameter<vector<double> >  ("ptHatLumi",vLumi);
  mPtHatBnd  = cfg.getUntrackedParameter<vector<double> >  ("ptHatBnd",vBnd);
  //------------- Compatibility check ------------------------------------------
  if (mRank != int(mMinPt.size())) {
    throw cms::Exception("DijetSearchHistos: ")<<"the min number of jets ("<<mRank<<") is not equal to the number of pt thresholds ("<<mMinPt.size()<<")";  
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void MultijetSearchGENHistos::beginJob() 
{
  mInf = TFile::Open(mFileName.c_str());
  mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());
  mTree = (TTree*)mDir->Get(mTreeName.c_str());
  mEvent = new QCDEvent();
  TBranch *branch = mTree->GetBranch("events");
  branch->SetAddress(&mEvent);
  char name[1000];
  //--------- book histos ---------------------------------
  TFileDirectory mPFDir   = fs->mkdir(mDirName+"pf");
  mPtHat      = mPFDir.make<TH1F>("PtHat","PtHat",700,0,3500);
  mPtHatAll   = mPFDir.make<TH1F>("PtHatAll","PtHatAll",700,0,3500);
  mHT         = mPFDir.make<TH1F>("HT","HT",200,0,5000);
  mHTAll      = mPFDir.make<TH1F>("HTAll","HTAll",200,0,5000);
  mHT4J       = mPFDir.make<TH1F>("HT4J","HT4J",200,0,5000);
  mHT4Jcut    = mPFDir.make<TH1F>("HT4Jcut","HT4Jcut",200,0,5000);
  mMAllJ      = mPFDir.make<TH1F>("MAllJ","MAllJ",140,0,7000);
  mM2J        = mPFDir.make<TH1F>("M2J","M2J",700,0,7000);
  mM4J        = mPFDir.make<TH1F>("M4J","M4J",700,0,7000);
  mM4Jcut     = mPFDir.make<TH1F>("M4Jcut","M4Jcut",700,0,7000);
  mDR         = mPFDir.make<TH1F>("DR","DR",100,0,10);
  mPtRatio    = mPFDir.make<TH1F>("PtRatio","PtRatio",100,0,1);
  mJetMulti   = mPFDir.make<TH1F>("JetMulti","JetMulti",20,0,20);
  mM4JvsHT4J  = mPFDir.make<TH2F>("M4JvsHT4J","M4JvsHT4J",200,0,5000,700,0,7000);
  mPtHat->Sumw2();
  mPtHatAll->Sumw2();
  mHT->Sumw2();
  mHT4J->Sumw2();
  mHT4Jcut->Sumw2();
  mHTAll->Sumw2();
  mMAllJ->Sumw2();
  mM4J->Sumw2();
  mM4Jcut->Sumw2();
  mM2J->Sumw2();
  mDR->Sumw2();
  mPtRatio->Sumw2();
  mJetMulti->Sumw2(); 
  mM4JvsHT4J->Sumw2();
  for(int j=0;j<8;j++) {
    sprintf(name,"JetPt%d",j);
    mPt[j] = mPFDir.make<TH1F>(name,name,350,0,3500);
    mPt[j]->Sumw2();
    sprintf(name,"JetEta%d",j);
    mEta[j] = mPFDir.make<TH1F>(name,name,24,-3,3);
    mEta[j]->Sumw2();
    sprintf(name,"JetPhi%d",j); 
    mPhi[j] = mPFDir.make<TH1F>(name,name,25,-3.15,3.15);
    mPhi[j]->Sumw2();
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void MultijetSearchGENHistos::endJob() 
{
  mInf->Close();
}
//////////////////////////////////////////////////////////////////////////////////////////
void MultijetSearchGENHistos::analyze(edm::Event const& evt, edm::EventSetup const& iSetup) 
{ 
  unsigned NEntries = mTree->GetEntries();
  cout<<"File: "<<mFileName<<endl;
  cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
  int decade = 0;
  int counter_kin(0),counter_ht(0);
  //---------- loop over the events ----------------------
  unsigned NN = NEntries;
  if (mNEvents > -1)
    NN = (unsigned)mNEvents;
  for(unsigned i=0;i<NN;i++) {
    double progress = 10.0*i/(1.0*NN);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;          
    mTree->GetEntry(i);
    //----------- mc weight -------------------------------- 
    double pthat = mEvent->evtHdr().pthat();
    double wt = mEvent->evtHdr().weight();
    if (mPtHatBnd.size() > 1) {
      int ipthat = getBin(pthat,mPtHatBnd);
      if (ipthat > -1) {
        double Leff = mPtHatLumi[ipthat];
        wt = 1./Leff; 
      }
      else {
        wt = 0.0;
      }
    }
    mPtHatAll->Fill(pthat,wt);
    //-------- cut flow ------------------------------------
    if (int(mEvent->nGenJets()) > mRank-1) {
      bool cutPt(true),cutEta(true),cutHT(true);
      double HT(0),HTAll(0);
      LorentzVector P4[8],P4AllJ(0,0,0,0);
      for (int j=0;j<mRank;j++) { 
        cutPt  *= (mEvent->genjet(j).pt() >= mMinPt[j]);
        cutEta *= (fabs(mEvent->genjet(j).eta()) <= mMaxEta);
        HT += mEvent->genjet(j).pt();
        P4[j] = mEvent->genjet(j);
        P4AllJ += P4[j];
      }
      cutHT = (HT >= mMinHT);
      int njets(mRank);
      HTAll = HT;
      for(unsigned k=mRank;k<mEvent->nGenJets();k++) {
        if (mEvent->genjet(k).pt() >= mMinPt[k] && fabs(mEvent->genjet(k).eta()) <= mMaxEta) {
          njets++;
          HTAll += mEvent->genjet(k).pt();
        }
      }
      if (cutPt && cutEta) {
        counter_kin++;
        if (cutHT) {
          counter_ht++;
          mMAllJ->Fill(P4AllJ.mass(),wt);
          mHT->Fill(HT,wt); 
          mHTAll->Fill(HTAll,wt);
          mPtHat->Fill(pthat);
          mJetMulti->Fill(njets,wt);
          //----- 4J combinatorics --------
          for(int j1=0;j1<mRank-3;j1++) {
            for(int j2=j1+1;j2<mRank-2;j2++) {
              for(int j3=j2+1;j3<mRank-1;j3++) {
                for(int j4=j3+1;j4<mRank;j4++) {
                  double m = (P4[j1]+P4[j2]+P4[j3]+P4[j4]).mass();
                  mM4J->Fill(m,wt);
                  double HT4J = P4[j1].pt()+P4[j2].pt()+P4[j3].pt()+P4[j4].pt();
                  mHT4J->Fill(HT4J,wt);
                  mM4JvsHT4J->Fill(HT4J,m,wt);
                  if (m > mScale4J*(HT4J-mOffset4J)) {
                    mM4Jcut->Fill(m,wt);
                    mHT4Jcut->Fill(HT4J,wt);
                  }
                }
              }
            }
          }
          //----- 2J combinatorics -------- 
          for(int j1=0;j1<mRank-1;j1++) {
            for(int j2=j1+1;j2<mRank;j2++) {
              double m2j = (P4[j1]+P4[j2]).mass();
              double dR  = sqrt(pow(P4[j1].eta()-P4[j2].eta(),2)+pow(P4[j1].phi()-P4[j2].phi(),2));     
              mM2J->Fill(m2j,wt);
              mDR->Fill(dR,wt); 
            }
          }
          mPtRatio->Fill(mEvent->genjet(mRank-1).pt()/mEvent->genjet(0).pt(),wt);
          for(int j=0;j<mRank;j++) {
            mPt[j]->Fill(mEvent->genjet(j).pt(),wt);
            mPhi[j]->Fill(mEvent->genjet(j).phi(),wt);
            mEta[j]->Fill(mEvent->genjet(j).eta(),wt);
          }// genjet loop
        }// ht cut
      }// kinematic cuts
    }// at least 8 pfjets
  }// tree loop
  cout<<"Events passing kinematic cuts:    "<<counter_kin<<endl;
  cout<<"Events passing the HT cut:        "<<counter_ht<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
MultijetSearchGENHistos::~MultijetSearchGENHistos() 
{
}
//////////////////////////////////////////////////////////////////////////////////////////
int MultijetSearchGENHistos::getBin(double x, const std::vector<double>& boundaries)
{
  int i;
  int n = boundaries.size()-1;
  if (x<boundaries[0] || x>=boundaries[n])
    return -1;
  for(i=0;i<n;i++) {
    if (x>=boundaries[i] && x<boundaries[i+1])
      return i;
   }
  return 0;
}
DEFINE_FWK_MODULE(MultijetSearchGENHistos);
