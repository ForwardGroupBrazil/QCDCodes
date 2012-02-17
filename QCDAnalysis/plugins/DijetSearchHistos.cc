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

#include "KKousour/QCDAnalysis/plugins/DijetSearchHistos.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

DijetSearchHistos::DijetSearchHistos(edm::ParameterSet const& cfg) 
{
  mMassBND   = cfg.getParameter<std::vector<double> > ("massBnd");
  mMinMass   = cfg.getParameter<double> ("minMass");
  mMinPt     = cfg.getParameter<double> ("minPt");
  mMaxEta    = cfg.getParameter<double> ("maxEta");
  mMaxDeta   = cfg.getParameter<double> ("maxDeta");
  mChiIN     = cfg.getParameter<double> ("chiIN");
  mChiOUT     = cfg.getParameter<double> ("chiOUT");
  mFileName  = cfg.getParameter<std::string> ("filename");
  mTreeName  = cfg.getParameter<std::string> ("treename");
  mDirName   = cfg.getParameter<std::string> ("dirname");
  mTriggers  = cfg.getParameter<std::vector<std::string> > ("triggers");
  mIsMC      = cfg.getParameter<bool> ("isMC");
  mNEvents   = cfg.getParameter<int>  ("nEvents"); 
  mJetID     = cfg.getParameter<int>  ("jetID");
  mHCALNoise = cfg.getParameter<int>  ("hcalNoiseFilter");
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetSearchHistos::beginJob() 
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
  TFileDirectory mCaloDir = fs->mkdir(mDirName+"calo");
  double auxMass[500];
  for(unsigned im=0;im<mMassBND.size();im++) {
    auxMass[im] = mMassBND[im];
  }
  sprintf(name,"METovSUMET");
  mPFMETovSUMET   = mPFDir.make<TH1F>(name,name,100,0,1.0001);
  mCaloMETovSUMET = mCaloDir.make<TH1F>(name,name,100,0,1.0001);
  
  sprintf(name,"JetPt");
  mPFPt = mPFDir.make<TH1F>(name,name,3500,0,3500);
  mCaloPt = mCaloDir.make<TH1F>(name,name,3500,0,3500);
  mPFPt->Sumw2();
  sprintf(name,"PtJJ");
  mPFPtJJ = mPFDir.make<TH1F>(name,name,3500,0,3500);
  mCaloPtJJ = mCaloDir.make<TH1F>(name,name,3500,0,3500);
  mPFPtJJ->Sumw2();
  mCaloPtJJ->Sumw2();
  sprintf(name,"JetEta");
  mPFEta = mPFDir.make<TH1F>(name,name,300,-3,3);
  mCaloEta = mCaloDir.make<TH1F>(name,name,300,-3,3);
  mPFEta->Sumw2();
  mCaloEta->Sumw2();
  sprintf(name,"EtaBoost");
  mPFEtaBoost = mPFDir.make<TH1F>(name,name,300,-3,3);
  mCaloEtaBoost = mCaloDir.make<TH1F>(name,name,300,-3,3);
  mPFEtaBoost->Sumw2();
  mCaloEtaBoost->Sumw2();
  sprintf(name,"Deta");
  mPFDeta = mPFDir.make<TH1F>(name,name,300,-3,3);
  mCaloDeta = mCaloDir.make<TH1F>(name,name,300,-3,3);
  mPFDeta->Sumw2();
  mCaloDeta->Sumw2();
  sprintf(name,"Dphi");
  mPFDphi = mPFDir.make<TH1F>(name,name,200,0,3.15);
  mCaloDphi = mCaloDir.make<TH1F>(name,name,200,0,3.15);
  mPFDphi->Sumw2();
  mCaloDphi->Sumw2();
  sprintf(name,"Mass");
  mPFMass = mPFDir.make<TH1F>(name,name,mMassBND.size()-1,auxMass);
  mCaloMass = mCaloDir.make<TH1F>(name,name,mMassBND.size()-1,auxMass);
  mPFMass->Sumw2();
  mCaloMass->Sumw2();
  sprintf(name,"MassIN");
  mPFMassIN = mPFDir.make<TH1F>(name,name,mMassBND.size()-1,auxMass);
  mCaloMassIN = mCaloDir.make<TH1F>(name,name,mMassBND.size()-1,auxMass);
  mPFMassIN->Sumw2();
  mCaloMassIN->Sumw2();
  sprintf(name,"MassOUT");
  mPFMassOUT = mPFDir.make<TH1F>(name,name,mMassBND.size()-1,auxMass);
  mCaloMassOUT = mCaloDir.make<TH1F>(name,name,mMassBND.size()-1,auxMass);
  mPFMassOUT->Sumw2();
  mCaloMassOUT->Sumw2();
  sprintf(name,"Chi");
  mPFChi = mPFDir.make<TH1F>(name,name,100,0,20);
  mCaloChi = mCaloDir.make<TH1F>(name,name,100,0,20);
  mPFChi->Sumw2();
  mCaloChi->Sumw2();
  sprintf(name,"ChiIN");
  mPFChiIN = mPFDir.make<TH1F>(name,name,100,0,20);
  mCaloChiIN = mCaloDir.make<TH1F>(name,name,100,0,20);
  mPFChiIN->Sumw2();
  mCaloChiIN->Sumw2();
  sprintf(name,"ChiOUT");
  mPFChiOUT = mPFDir.make<TH1F>(name,name,100,0,20);
  mCaloChiOUT = mCaloDir.make<TH1F>(name,name,100,0,20);
  mPFChiOUT->Sumw2();
  mCaloChiOUT->Sumw2();

  //-------- jet properties ----------------------------
  sprintf(name,"CHF");
  mCHF = mPFDir.make<TH1F>(name,name,100,0,1.001);
  mCHF->Sumw2();  
  sprintf(name,"NHF");
  mNHF = mPFDir.make<TH1F>(name,name,100,0,1.001);
  mNHF->Sumw2();
  sprintf(name,"PHF");
  mPHF = mPFDir.make<TH1F>(name,name,100,0,1.001);
  mPHF->Sumw2();

  sprintf(name,"CHFVsRun"); 
  mCHFVsRun = mPFDir.make<TProfile>(name,name,3,0,3,0,1.1);
  mCHFVsRun->SetBit(TH1::kCanRebin);
  sprintf(name,"NHFVsRun"); 
  mNHFVsRun = mPFDir.make<TProfile>(name,name,3,0,3,0,1.1);
  mNHFVsRun->SetBit(TH1::kCanRebin);
  sprintf(name,"PHFVsRun"); 
  mPHFVsRun = mPFDir.make<TProfile>(name,name,3,0,3,0,1.1);
  mPHFVsRun->SetBit(TH1::kCanRebin);
      
  sprintf(name,"EMF");
  mEMF = mCaloDir.make<TH1F>(name,name,100,0,1.001);
  mEMF->Sumw2();
  sprintf(name,"N90hits");
  mN90hits = mCaloDir.make<TH1F>(name,name,200,0,200);
  mN90hits->Sumw2();
  sprintf(name,"fHPD");
  mfHPD = mCaloDir.make<TH1F>(name,name,100,0,1.001);
  mfHPD->Sumw2();
  sprintf(name,"NTrkCalo");
  mNTrkCalo = mCaloDir.make<TH1F>(name,name,200,0,200);
  mNTrkCalo->Sumw2();
  sprintf(name,"NTrkVtx");
  mNTrkVtx = mCaloDir.make<TH1F>(name,name,200,0,200);
  mNTrkVtx->Sumw2();
      
  sprintf(name,"EMFVsRun"); 
  mEMFVsRun = mCaloDir.make<TProfile>(name,name,3,0,3,0,1.1);
  mEMFVsRun->SetBit(TH1::kCanRebin);
  sprintf(name,"NTrkCaloVsRun");
  mNTrkCaloVsRun = mCaloDir.make<TProfile>(name,name,3,0,3,0,1.1);
  mNTrkCaloVsRun->SetBit(TH1::kCanRebin);
  sprintf(name,"NTrkVtxVsRun");
  mNTrkVtxVsRun = mCaloDir.make<TProfile>(name,name,3,0,3,0,1.1);
  mNTrkVtxVsRun->SetBit(TH1::kCanRebin);
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetSearchHistos::endJob() 
{
  mInf->Close();
}
//////////////////////////////////////////////////////////////////////////////////////////
int DijetSearchHistos::findRun(int x, const std::vector<int>& runs)
{
  int result(-1);
  for(unsigned i=0;i<runs.size();i++)
    if (x == runs[i])
      return i;
  return result;
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetSearchHistos::analyze(edm::Event const& evt, edm::EventSetup const& iSetup) 
{ 
  unsigned NEntries = mTree->GetEntries();
  cout<<"File: "<<mFileName<<endl;
  cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
  int decade = 0;
  int counter_hlt(0),counter_pv(0),counter_hcal(0);
  int pf_counter_all(0),calo_counter_all(0);
  int Ntrig = (int)mTriggers.size();
  vector<int> Runs;
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
    bool hltPass(false);
    for(int itrig=0;itrig<TMath::Max(1,Ntrig);itrig++) {
      if (Ntrig == 0) {
        hltPass = true;
      }
      else {
        int ihlt = mTrigIndex[itrig];
        if (mEvent->fired(ihlt) > 0) {
          hltPass = true;
          continue;
        }
      }
    }
    //-------- cut flow ------------------------------------
    if (hltPass) {
      counter_hlt++;
      if (mEvent->evtHdr().isPVgood() == 1) {
        counter_pv++;
        bool cut_hcalNoise(true);
        cut_hcalNoise = mEvent->evtHdr().hcalNoise();
        if (cut_hcalNoise) {
          counter_hcal++;  
          if (mEvent->nPFJets() > 1) {
            double deta = mEvent->pfjet(0).eta()-mEvent->pfjet(1).eta();
            double etaBoost = 0.5*(mEvent->pfjet(0).eta()+mEvent->pfjet(1).eta());
            double ptJJ = (mEvent->pfjet(0).p4()+mEvent->pfjet(1).p4()).pt();  
            double dphi = mEvent->pfjet(0).phi()-mEvent->pfjet(1).phi();
            double mass = mEvent->pfmjjcor(0);
            double cosThetaStar = tanh(0.5*deta);
            double chi = (1+fabs(cosThetaStar))/(1-fabs(cosThetaStar));
            bool cutPt   = (mEvent->pfjet(0).ptCor() >= mMinPt && mEvent->pfjet(1).ptCor() >= mMinPt);
            bool cutEta  = (fabs(mEvent->pfjet(0).eta()) <= mMaxEta && fabs(mEvent->pfjet(1).eta()) <= mMaxEta);
            bool cutDeta = (fabs(deta) <= mMaxDeta);
            bool cutMass = (mass >= mMinMass);
            bool cutID(true);
            bool cutIN = (chi <= mChiIN);
            bool cutOUT = (chi > mChiIN && chi <=mChiOUT);
            if (mJetID == 1)
              cutID = (mEvent->pfjet(0).looseID() && mEvent->pfjet(1).looseID());
            if (mJetID == 2)
              cutID = (mEvent->pfjet(0).tightID() && mEvent->pfjet(1).tightID());  
            if (cutPt && cutEta && cutID && cutMass) {
              mPFPtJJ->Fill(ptJJ,wt);
              mPFEtaBoost->Fill(etaBoost,wt);
              mPFChi->Fill(chi,wt);
              if (cutIN) {
                mPFMassIN->Fill(mass,wt);
                mPFChiIN->Fill(chi,wt);
              }
              if (cutOUT) {
                mPFMassOUT->Fill(mass,wt);
                mPFChiOUT->Fill(chi,wt);
              }
              if (cutDeta) {
                pf_counter_all++;
                mPFMass->Fill(mass,wt);
                mPFDeta->Fill(deta,wt);
                mPFDphi->Fill(fabs(dphi),wt);
                mPFMETovSUMET->Fill(mEvent->pfmet().met_o_sumet(),wt);
                for(unsigned j=0;j<2;j++) {
                  mPFPt->Fill((mEvent->pfjet(j)).ptCor(),wt);
                  mPFEta->Fill((mEvent->pfjet(j)).eta(),wt); 
                  mCHF->Fill((mEvent->pfjet(j)).chf(),wt);
                  mNHF->Fill((mEvent->pfjet(j)).nhf(),wt);
                  mPHF->Fill((mEvent->pfjet(j)).phf(),wt); 
                  sprintf(name,"%d",mEvent->evtHdr().runNo());	
                  mCHFVsRun->Fill(name,(mEvent->pfjet(j)).chf(),1);
                  mNHFVsRun->Fill(name,(mEvent->pfjet(j)).nhf(),1); 
                  mPHFVsRun->Fill(name,(mEvent->pfjet(j)).phf(),1);
                }// pfjet loop
              }// cut deta
            }// cuts
          }// at least two pfjets
          if (mEvent->nCaloJets() > 1) {
            double deta = mEvent->calojet(0).eta()-mEvent->calojet(1).eta();
            double dphi = mEvent->calojet(0).phi()-mEvent->calojet(1).phi();
            double etaBoost = 0.5*(mEvent->calojet(0).eta()+mEvent->calojet(1).eta());
            double ptJJ = (mEvent->calojet(0).p4()+mEvent->calojet(1).p4()).pt();
            double mass = mEvent->calomjjcor(0);
            double cosThetaStar = tanh(0.5*deta);
            double chi = (1+fabs(cosThetaStar))/(1-fabs(cosThetaStar));
            bool cutPt   = (mEvent->calojet(0).ptCor() >= mMinPt && mEvent->calojet(1).ptCor() >= mMinPt);
            bool cutEta  = (fabs(mEvent->calojet(0).eta()) <= mMaxEta && fabs(mEvent->calojet(1).eta()) <= mMaxEta); 
            bool cutDeta = (fabs(deta) <= mMaxDeta);
            bool cutMass = (mass >= mMinMass);
            bool cutID(true);
            bool cutIN = (chi <= mChiIN);
            bool cutOUT = (chi > mChiIN && chi <= mChiOUT);
            if (mJetID == 1)
              cutID = (mEvent->calojet(0).looseID() && mEvent->calojet(1).looseID());
            if (mJetID == 2)
              cutID = (mEvent->calojet(0).tightID() && mEvent->calojet(1).tightID());
            if (cutPt && cutEta && cutID && cutMass) {
              mCaloPtJJ->Fill(ptJJ,wt);
              mCaloEtaBoost->Fill(etaBoost,wt);
              mCaloChi->Fill(chi,wt);
              if (cutIN) {
                mCaloMassIN->Fill(mass,wt);
                mCaloChiIN->Fill(chi,wt);
              }
              if (cutOUT) {
                mCaloMassOUT->Fill(mass,wt);
                mCaloChiOUT->Fill(chi,wt);
              }  
              if (cutDeta) {
                calo_counter_all++;
                mCaloMass->Fill(mass,wt);
                mCaloDeta->Fill(deta,wt);
                mCaloDphi->Fill(fabs(dphi),wt);
                mCaloMETovSUMET->Fill(mEvent->calomet().met_o_sumet(),wt);
                for(unsigned j=0;j<2;j++) {
                  mCaloPt->Fill((mEvent->calojet(j)).ptCor(),wt);
                  mCaloEta->Fill((mEvent->calojet(j)).eta(),wt);
                  mEMF->Fill((mEvent->calojet(j)).emf(),wt);
                  mN90hits->Fill((mEvent->calojet(j)).n90hits(),wt);
                  mfHPD->Fill((mEvent->calojet(j)).fHPD(),wt);
                  mNTrkCalo->Fill((mEvent->calojet(j)).nTrkCalo(),wt);
                  mNTrkVtx->Fill((mEvent->calojet(j)).nTrkVtx(),wt);
                  sprintf(name,"%d",mEvent->evtHdr().runNo());
                  mEMFVsRun->Fill(name,(mEvent->calojet(j)).emf(),1);
                  mNTrkCaloVsRun->Fill(name,(mEvent->calojet(j)).nTrkCalo(),1);
                  mNTrkVtxVsRun->Fill(name,(mEvent->calojet(j)).nTrkVtx(),1);  
                }// calojet loop
              }// cut deta
            }// cuts 
          }// at least two calojets 
        }// if pv
      }// if hlt
    }// trigger loop
  }// tree loop
  if (Ntrig > 0) {
    cout<<"Events passing the trigger:    "<<counter_hlt<<endl;
  } 
  cout<<"Events passing the PV:         "<<counter_pv<<endl;
  cout<<"Events passing the HCAL NOISE: "<<counter_hcal<<endl;  
  cout<<"Events passing all PF cuts:    "<<pf_counter_all<<endl;
  cout<<"Events passing all CALO cuts:  "<<calo_counter_all<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
DijetSearchHistos::~DijetSearchHistos() 
{
}

DEFINE_FWK_MODULE(DijetSearchHistos);
