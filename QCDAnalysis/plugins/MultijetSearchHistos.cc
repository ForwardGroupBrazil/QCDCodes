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

#include "KKousour/QCDAnalysis/plugins/MultijetSearchHistos.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

MultijetSearchHistos::MultijetSearchHistos(edm::ParameterSet const& cfg) 
{
  mFileName  = cfg.getParameter<std::string>               ("filename");
  mTreeName  = cfg.getParameter<std::string>               ("treename");
  mDirName   = cfg.getParameter<std::string>               ("dirname");
  mRank      = cfg.getParameter<int>                       ("rank");
  mNEvents   = cfg.getParameter<int>                       ("nEvents"); 
  mJetID     = cfg.getParameter<int>                       ("jetID");
  mHCALNoise = cfg.getParameter<int>                       ("hcalNoiseFilter");
  mMaxEta    = cfg.getParameter<double>                    ("maxEta");
  mMinPt     = cfg.getParameter<std::vector<double> >      ("minPt");
  mTriggers  = cfg.getParameter<std::vector<std::string> > ("triggers");
  //------------- MC -----------------------------------------------------------
  vector<double> vLumi(0),vBnd(0);
  mIsMC      = cfg.getUntrackedParameter<bool>             ("isMC",false);
  mPtHatLumi = cfg.getUntrackedParameter<vector<double> >  ("ptHatLumi",vLumi);
  mPtHatBnd  = cfg.getUntrackedParameter<vector<double> >  ("ptHatBnd",vBnd);
  //------------- Compatibility check ------------------------------------------
  if (mRank != int(mMinPt.size())) {
    throw cms::Exception("DijetSearchHistos: ")<<"the min number of jets ("<<mRank<<") is not equal to the number of pt thresholds ("<<mMinPt.size()<<")";  
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void MultijetSearchHistos::beginJob() 
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
        throw cms::Exception("DijetSearchHistos: ")<<"The requested trigger ("<<mTriggers[itrig]<<") is not found ";
      }
    }
    for(unsigned itrig=0;itrig<mTriggers.size();itrig++)
    cout<<mTriggers[itrig]<<" --> "<<mTrigIndex[itrig]<<endl;
  } 
  //--------- book histos ---------------------------------
  TFileDirectory mPFDir   = fs->mkdir(mDirName+"pf");
  mMETovSUMET = mPFDir.make<TH1F>("METovSUMET","METovSUMET",100,0,1.0001);
  mPtHat      = mPFDir.make<TH1F>("PtHat","PtHat",700,0,3500);
  mPtHatAll   = mPFDir.make<TH1F>("PtHatAll","PtHatAll",700,0,3500);
  mHT         = mPFDir.make<TH1F>("HT","HT",200,0,5000);
  mHTAll      = mPFDir.make<TH1F>("HTAll","HTAll",200,0,5000);
  mMAllJ      = mPFDir.make<TH1F>("MAllJ","MAllJ",140,0,7000);
  mM2J        = mPFDir.make<TH1F>("M2J","M2J",700,0,7000);
  mM4J        = mPFDir.make<TH1F>("M4J","M4J",700,0,7000);
  mDR         = mPFDir.make<TH1F>("DR","DR",100,0,10);
  mJetMulti   = mPFDir.make<TH1F>("JetMulti","JetMulti",20,0,20);
  mMETovSUMET->Sumw2();
  mPtHat->Sumw2();
  mPtHatAll->Sumw2();
  mHT->Sumw2();
  mHTAll->Sumw2();
  mMAllJ->Sumw2();
  mM4J->Sumw2();
  mM2J->Sumw2();
  mDR->Sumw2();
  mJetMulti->Sumw2(); 
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
    //-------- jet properties ----------------------------
    sprintf(name,"CHF%d",j);
    mCHF[j] = mPFDir.make<TH1F>(name,name,100,0,1.001);
    mCHF[j]->Sumw2();  
    sprintf(name,"NHF%d",j);
    mNHF[j] = mPFDir.make<TH1F>(name,name,100,0,1.001);
    mNHF[j]->Sumw2();
    sprintf(name,"PHF%d",j);
    mPHF[j] = mPFDir.make<TH1F>(name,name,100,0,1.001);
    mPHF[j]->Sumw2();
    sprintf(name,"ELF%d",j);
    mELF[j] = mPFDir.make<TH1F>(name,name,100,0,1.001);
    mELF[j]->Sumw2();
    sprintf(name,"MUF%d",j);
    mMUF[j] = mPFDir.make<TH1F>(name,name,100,0,1.001);
    mMUF[j]->Sumw2();
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void MultijetSearchHistos::endJob() 
{
  mInf->Close();
}
//////////////////////////////////////////////////////////////////////////////////////////
void MultijetSearchHistos::analyze(edm::Event const& evt, edm::EventSetup const& iSetup) 
{ 
  unsigned NEntries = mTree->GetEntries();
  cout<<"File: "<<mFileName<<endl;
  cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
  int decade = 0;
  int counter_hlt(0),counter_pv(0),counter_hcal(0),counter_kin(0),counter_id(0);
  int Ntrig = (int)mTriggers.size();
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
    double wt(1.0),pthat(0.0); 
    if (mIsMC) {
      pthat = mEvent->evtHdr().pthat();
      wt = mEvent->evtHdr().weight();
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
    }
    //---------- loop over the triggers ----------------------
    bool hltPass(false);
    int prescale(1);
    for(int itrig=0;itrig<TMath::Max(1,Ntrig);itrig++) {
      if (Ntrig == 0) {
        hltPass = true;
      }
      else {
        int ihlt = mTrigIndex[itrig];
        if (mEvent->fired(ihlt) > 0) {
          hltPass = true;
          prescale = mEvent->preL1(ihlt) * mEvent->preHLT(ihlt);
          continue;
        }
      }
    }
    wt *= prescale;
    //-------- cut flow ------------------------------------
    if (hltPass) {
      counter_hlt++;
      if (mEvent->evtHdr().isPVgood() == 1) {
        counter_pv++;
        bool cut_hcalNoise(true);
        if (mHCALNoise == 1)
          cut_hcalNoise = mEvent->evtHdr().looseHCALNoise();
        if (mHCALNoise == 2)
          cut_hcalNoise = mEvent->evtHdr().tightHCALNoise(); 
        if (cut_hcalNoise) {
          counter_hcal++;  
          if (int(mEvent->nPFJets()) > mRank-1) {
            bool cutPt(true),cutEta(true),cutID(true);
            double HT(0),HTAll(0);
            LorentzVector P4[8],P4AllJ(0,0,0,0);
            for (int j=0;j<mRank;j++) { 
              cutPt  *= (mEvent->pfjet(j).ptCor() >= mMinPt[j]);
              cutEta *= (fabs(mEvent->pfjet(j).eta()) <= mMaxEta);
              if (mJetID == 1)//---loose
                cutID *= (mEvent->pfjet(j).looseID());
              if (mJetID == 2)//---tight
                cutID *= (mEvent->pfjet(j).tightID());  
              if (mJetID == 3)//---super tight
                cutID *= (mEvent->pfjet(j).tightID() && mEvent->pfjet(j).muf() < 0.9 && mEvent->pfjet(j).elf() < 0.9);
              HT += mEvent->pfjet(j).ptCor();
              P4[j] = mEvent->pfjet(j).cor() * mEvent->pfjet(j).p4();
              P4AllJ += P4[j];
            }
            int njets(mRank);
            HTAll = HT;
            for(unsigned k=mRank;k<mEvent->nPFJets();k++) {
              bool id(false);
              if (mJetID == 1)
                id = (mEvent->pfjet(k).looseID());
              if (mJetID == 2)
                id = (mEvent->pfjet(k).tightID());
              if (mJetID == 3)
                id = (mEvent->pfjet(k).tightID() && mEvent->pfjet(k).muf() < 0.9 && mEvent->pfjet(k).elf() < 0.9);
              if (mEvent->pfjet(k).ptCor() >= mMinPt[mRank-1] && fabs(mEvent->pfjet(k).eta()) <= mMaxEta && id) { 
                njets++;
                HTAll += mEvent->pfjet(k).ptCor();
              }
            } 
            if (cutPt && cutEta) {
              counter_kin++;
              if (cutID) {
                counter_id++;
                mMAllJ->Fill(P4AllJ.mass(),wt);
                mMETovSUMET->Fill(mEvent->pfmet().met_o_sumet(),wt);
                mHT->Fill(HT,wt); 
                mHTAll->Fill(HTAll,wt);
                mPtHat->Fill(pthat);
                mJetMulti->Fill(njets,wt);
                //----- 4J combinatorics --------
                for(int j1=0;j1<mRank-3;j1++) {
                  for(int j2=j1+1;j2<mRank-2;j2++) {
                    for(int j3=j2+1;j3<mRank-1;j3++) {
                      for(int j4=j3+1;j4<mRank;j4++) {
                        mM4J->Fill((P4[j1]+P4[j2]+P4[j3]+P4[j4]).mass(),wt);
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
                for(int j=0;j<mRank;j++) {
                  mPt[j]->Fill((mEvent->pfjet(j)).ptCor(),wt);
                  mPhi[j]->Fill((mEvent->pfjet(j)).phi(),wt);
                  mEta[j]->Fill((mEvent->pfjet(j)).eta(),wt);
                  mCHF[j]->Fill((mEvent->pfjet(j)).chf(),wt);
                  mNHF[j]->Fill((mEvent->pfjet(j)).nhf(),wt);
                  mPHF[j]->Fill((mEvent->pfjet(j)).phf(),wt);
                  mELF[j]->Fill((mEvent->pfjet(j)).elf(),wt);
                  mMUF[j]->Fill((mEvent->pfjet(j)).muf(),wt);
                }// pfjet loop
              }// id cuts
            }// kinematic cuts
          }// at least 8 pfjets
        }// if pv
      }// if hlt
    }// trigger loop
  }// tree loop
  if (Ntrig > 0) {
    cout<<"Events passing the trigger:       "<<counter_hlt<<endl;
  } 
  cout<<"Events passing the PV:            "<<counter_pv<<endl;
  cout<<"Events passing the HCAL NOISE:    "<<counter_hcal<<endl;  
  cout<<"Events passing kinematic cuts:    "<<counter_kin<<endl;
  cout<<"Events passing the jet id cuts:   "<<counter_id<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
MultijetSearchHistos::~MultijetSearchHistos() 
{
}
//////////////////////////////////////////////////////////////////////////////////////////
int MultijetSearchHistos::getBin(double x, const std::vector<double>& boundaries)
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
DEFINE_FWK_MODULE(MultijetSearchHistos);