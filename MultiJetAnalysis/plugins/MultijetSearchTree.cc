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

#include "KKousour/MultiJetAnalysis/plugins/MultijetSearchTree.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

MultijetSearchTree::MultijetSearchTree(edm::ParameterSet const& cfg) 
{
  mFileNames   = cfg.getParameter<std::vector<std::string> > ("filenames");
  mTreeName    = cfg.getParameter<std::string>               ("treename");
  mDirName     = cfg.getParameter<std::string>               ("dirname");
  mTriggers    = cfg.getParameter<std::vector<std::string> > ("triggers");
  mIsMC        = cfg.getUntrackedParameter<bool>             ("isMC",false);
  mIsPreScaled = cfg.getUntrackedParameter<bool>             ("isPreScaled",false);
  vector<double> vLumi(0),vBnd(0);
  mPtHatLumi   = cfg.getUntrackedParameter<vector<double> >  ("ptHatLumi",vLumi);
  mPtHatBnd    = cfg.getUntrackedParameter<vector<double> >  ("ptHatBnd",vBnd); 
  mNEvents     = cfg.getParameter<int>                       ("nEvents"); 
  mJetID       = cfg.getParameter<int>                       ("jetID");
  mHCALNoise   = cfg.getParameter<int>                       ("hcalNoiseFilter");
  mEtaMax      = cfg.getParameter<double>                    ("etaMAX");
  mPtMin       = cfg.getParameter<double>                    ("ptMIN");
  mPUTagMin    = cfg.getParameter<double>                    ("puTagMIN");
}
//////////////////////////////////////////////////////////////////////////////////////////
void MultijetSearchTree::beginJob() 
{
  mInf = TFile::Open(mFileNames[0].c_str());
  mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());
  mEvent = new QCDEvent();
  mOutTree = fs->make<TTree>("tr","tr");
  mOutTree->Branch("runNo"       ,&mRun         ,"mRun/I");
  mOutTree->Branch("evtNo"       ,&mEvt         ,"mEvt/I");
  mOutTree->Branch("npv"         ,&mNPV         ,"mNPV/I");
  mOutTree->Branch("m4jIndex"    ,&mM4JIndex    ,"mM4JIndex[2]/I");
  mOutTree->Branch("m2jIndex"    ,&mM2JIndex    ,"mM2JIndex[2][2]/I");
  mOutTree->Branch("dphi4j"      ,&mDphi4J      ,"mDphi4J/F");
  mOutTree->Branch("ht"          ,&mHT          ,"mHT/F");
  mOutTree->Branch("m8j"         ,&mM8J         ,"mM8J/F");
  mOutTree->Branch("weight"      ,&mWeight      ,"mWeight/F");
  mOutTree->Branch("cosThetaStar",&mCosThetaStar,"mCosThetaStar/F");
  mOutTree->Branch("m4jBalance"  ,&mM4JBalance  ,"mM4JBalance/F");
  mOutTree->Branch("m2jBalance"  ,&mM2JBalance  ,"mM2JBalance[2]/F"); 
  mOutTree->Branch("m4j"         ,&mM4J         ,"mM4J[2]/F");
  mOutTree->Branch("ht4j"        ,&mHT4J        ,"mHT4J[2]/F");
  mOutTree->Branch("m2j"         ,&mM2J         ,"mM2J[2][2]/F");
  mOutTree->Branch("m2jAll"      ,&mM2JAll      ,"mM2JAll[28]/F");
  mOutTree->Branch("dR2jAll"     ,&mDR2JAll     ,"mDR2JAll[28]/F");
  mOutTree->Branch("m4jAll"      ,&mM4JAll      ,"mM4JAll[70]/F");
  mOutTree->Branch("ht4jAll"     ,&mHT4JAll     ,"mHT4JAll[70]/F");
  mOutTree->Branch("pt"          ,&mPt          ,"mPt[8]/F");
  mOutTree->Branch("beta"        ,&mBeta        ,"mBeta[8]/F");
  mOutTree->Branch("betaStar"    ,&mBetaStar    ,"mBetaStar[8]/F");
  mOutTree->Branch("puTag"       ,&mPUTag       ,"mPUTag[8]/F");
  mOutTree->Branch("eta"         ,&mEta         ,"mEta[8]/F");
  mOutTree->Branch("phi"         ,&mPhi         ,"mPhi[8]/F");
  mOutTree->Branch("mass"        ,&mMass        ,"mMass[8]/F");
  mOutTree->Branch("chf"         ,&mCHF         ,"mCHF[8]/F");
  mOutTree->Branch("nhf"         ,&mNHF         ,"mNHF[8]/F");
  mOutTree->Branch("phf"         ,&mPHF         ,"mPHF[8]/F");
  mOutTree->Branch("muf"         ,&mMUF         ,"mMUF[8]/F");
  mOutTree->Branch("elf"         ,&mELF         ,"mELF[8]/F");
  if (mIsMC) {
    mOutTree->Branch("pthat",&mPtHat,"mPtHat/F");
  }
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
        throw cms::Exception("MultiSearchTree: ")<<"The requested trigger ("<<mTriggers[itrig]<<") is not found ";
      }
    }
    for(unsigned itrig=0;itrig<mTriggers.size();itrig++)
    cout<<mTriggers[itrig]<<" --> "<<mTrigIndex[itrig]<<endl;
  } 
}
//////////////////////////////////////////////////////////////////////////////////////////
void MultijetSearchTree::endJob() 
{

}
//////////////////////////////////////////////////////////////////////////////////////////
int MultijetSearchTree::getBin(double x, const std::vector<double>& boundaries)
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
void MultijetSearchTree::analyze(edm::Event const& evt, edm::EventSetup const& iSetup) 
{
  int counter_hlt(0),counter_pv(0),counter_hcal(0),counter_id(0);
  int Ntrig = (int)mTriggers.size();
  for(unsigned iFile=0;iFile<mFileNames.size();iFile++) {
    TFile *inf = TFile::Open(mFileNames[iFile].c_str());
    TTree *mTree = (TTree*)inf->Get((mDirName+"/"+mTreeName).c_str());
    unsigned long NEntries = mTree->GetEntries();
    cout<<"Adding file: "<<mFileNames[iFile]<<endl;
    cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
    TBranch *branch = mTree->GetBranch("events");
    branch->SetAddress(&mEvent); 
    int decade = 0;
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
      double wt(1.0);
      if (mIsMC) {
        mPtHat = mEvent->evtHdr().pthat();
        wt = mEvent->evtHdr().weight();
        if (mPtHatBnd.size() > 1) {
          int ipthat = getBin(mPtHat,mPtHatBnd);
          if (ipthat > -1) {
            double Leff = mPtHatLumi[ipthat];
            wt = 1./Leff; 
          }
          else {
            wt = 0.0;
          }
        } 
      }
      //---------- loop over the triggers ----------------------
      bool hltPass(false);
      double prescale(1.0);
      for(int itrig=0;itrig<TMath::Max(1,Ntrig);itrig++) {
        if (Ntrig == 0) {
          hltPass = true;
        }
        else {
          int ihlt = mTrigIndex[itrig];
          if (mEvent->fired(ihlt) > 0) {
            prescale = mEvent->preL1(ihlt) * mEvent->preHLT(ihlt);
            if (mIsPreScaled) {
              hltPass = true;
              continue;
            }
            else {
              if (prescale == 1) {
                hltPass = true; 
                continue;
              }
            }
          }
        }
      }
      wt = prescale;
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
            double HT(0),HT4J_1(0),HT4J_2(0),M4J(0),M4J_1(0),M4J_2(0),cosThetaStar(0);
            if (mEvent->nPFJets() > 7) {
              bool cutID(true);
              bool cut_eta(true);
              bool cut_pt = (mEvent->pfjet(7).ptCor() > mPtMin);
              bool cut_pu(true);
              LorentzVector P4[8],P48J(0,0,0,0);
              for (int j=0;j<8;j++) { 
                cut_pu *= ((1-mEvent->pfjet(j).betaStar()) > mPUTagMin);
                if (mJetID == 1)
                  cutID *= (mEvent->pfjet(j).looseID());
                if (mJetID == 2)
                  cutID *= (mEvent->pfjet(j).tightID());  
                HT += mEvent->pfjet(j).ptCor();
                P4[j] = mEvent->pfjet(j).cor() * mEvent->pfjet(j).p4();
                P48J += P4[j];
                cut_eta *= (fabs(mEvent->pfjet(j).eta()) < mEtaMax);
              }
              if (cutID && cut_pu && cut_eta && cut_pt) {
                mM8J    = P48J.mass();
                mHT     = HT;
                mRun    = mEvent->evtHdr().runNo();
                mEvt    = mEvent->evtHdr().event();
                mNPV    = mEvent->evtHdr().nVtxGood();
                mWeight = wt;
                //----- 4J combinatorics --------
                int jj(0);
                for(int j1=0;j1<5;j1++) {
                  for(int j2=j1+1;j2<6;j2++) {
                    for(int j3=j2+1;j3<7;j3++) {
                      for(int j4=j3+1;j4<8;j4++) {
                        mM4JAll[jj]  = (P4[j1]+P4[j2]+P4[j3]+P4[j4]).mass();
                        mHT4JAll[jj] = P4[j1].pt()+P4[j2].pt()+P4[j3].pt()+P4[j4].pt();
                        jj++;
                      }
                    }
                  }
                }
                //----- most probable 4J combinations -----------             
                vector<int> v1,v2,v4J[2]; 
                LorentzVector P44J_1,P44J_2;
                double d,dmin(10000);
                for(int j1=0;j1<1;j1++) {
                  for(int j2=j1+1;j2<6;j2++) {
                    for(int j3=j2+1;j3<7;j3++) {
                      for(int j4=j3+1;j4<8;j4++) {
                        v1.clear();
                        v2.clear();
                        v1.push_back(j1);
                        v1.push_back(j2);
                        v1.push_back(j3);
                        v1.push_back(j4);   
                        // find the complementary 4 jet indices
                        v2 = findIndices(v1,4);
                        P44J_1 = P4[v1[0]]+P4[v1[1]]+P4[v1[2]]+P4[v1[3]];
                        P44J_2 = P4[v2[0]]+P4[v2[1]]+P4[v2[2]]+P4[v2[3]];
                        M4J_1  = P44J_1.mass();
                        M4J_2  = P44J_2.mass();
                        HT4J_1 = P4[v1[0]].pt()+P4[v1[1]].pt()+P4[v1[2]].pt()+P4[v1[3]].pt();
                        HT4J_2 = P4[v2[0]].pt()+P4[v2[1]].pt()+P4[v2[2]].pt()+P4[v2[3]].pt();
                        cosThetaStar = tanh(0.5*(P44J_1.eta()-P44J_2.eta()));
                        d = fabs(M4J_1-M4J_2);
                        if (d < dmin) {
                          mM4J[0] = M4J_1;
                          mM4J[1] = M4J_2;
                          M4J = 0.5*(M4J_1+M4J_2);
                          mHT4J[0] = HT4J_1;
                          mHT4J[1] = HT4J_2;
                          mCosThetaStar = cosThetaStar;
                          mDphi4J = deltaPhi(P44J_1.phi(),P44J_2.phi());  
                          mM4JIndex[0] = (v1[0]+1)*1000+(v1[1]+1)*100+(v1[2]+1)*10+(v1[3]+1);
                          mM4JIndex[1] = (v2[0]+1)*1000+(v2[1]+1)*100+(v2[2]+1)*10+(v2[3]+1);
                          v4J[0] = v1;
                          v4J[1] = v2;
                          dmin = d;
                        } 
                      }
                    }
                  }
                }
                mM4JBalance   = 2*fabs(mM4J[0]-mM4J[1])/(mM4J[0]+mM4J[1]);
                //----- 2J combinatorics -------- 
                jj = 0;
                for(int j1=0;j1<8;j1++) {
                  for(int j2=j1+1;j2<8;j2++) {
                    mM2JAll[jj] = (P4[j1]+P4[j2]).mass();
                    mDR2JAll[jj] = deltaR(P4[j1],P4[j2]);
                    jj++;
                  }
                } 
                LorentzVector P42J[2];
                double M2J[2];
                for(int k=0;k<2;k++) {// loop over the quartrets
                  dmin = 10000;
                  for(int j1=0;j1<1;j1++) {
                    for(int j2=j1+1;j2<4;j2++) {
                      v1.clear();
                      v2.clear();
                      v1.push_back(j1);
                      v1.push_back(j2);
                      v2 = findIndices(v1,2);
                      P42J[0] = P4[v4J[k][v1[0]]]+P4[v4J[k][v1[1]]];
                      P42J[1] = P4[v4J[k][v2[0]]]+P4[v4J[k][v2[1]]];
                      M2J[0] = P42J[0].mass();
                      M2J[1] = P42J[1].mass();
                      d = fabs(M2J[0]-M2J[1]);
                      if (d < dmin) {
               //cout<<k<<" Found new min: "<<v4J[k][v1[0]]<<" "<<v4J[k][v1[1]]<<" "<<v4J[k][v2[0]]<<" "<<v4J[k][v2[1]]<<" "<<M2J[0]<<" "<<M2J[1]<<endl; 
                        mM2J[k][0] = M2J[0];
                        mM2J[k][1] = M2J[1];
                        mM2JBalance[k] = 2*fabs(mM2J[k][0]-mM2J[k][1])/(mM2J[k][0]+mM2J[k][1]);
                        mM2JIndex[k][0] = (v4J[k][v1[0]]+1)*10+(v4J[k][v1[1]]+1);
                        mM2JIndex[k][1] = (v4J[k][v2[0]]+1)*10+(v4J[k][v2[1]]+1);
                        dmin = d;
                      }
                    }    
                  }
                }
                /*
                cout<<"First 4j group: "<<mM4JIndex[0]<<endl;
                cout<<"2j doublets: "<<mM2JIndex[0][0]<<" "<<mM2JIndex[0][1]<<endl;
                cout<<"Second 4j group: "<<mM4JIndex[1]<<endl;
                cout<<"2j doublets: "<<mM2JIndex[1][0]<<" "<<mM2JIndex[1][1]<<endl;
                */
                for(int j=0;j<8;j++) {
                  mPt[j]  = mEvent->pfjet(j).ptCor();
                  if (!mIsMC) {
                    //--- to be fixed when beta will be available in the MC ntuples
                    mBeta[j]     = mEvent->pfjet(j).beta();
                    mBetaStar[j] = mEvent->pfjet(j).betaStar();
                    mPUTag[j]    = 1-mBetaStar[j];
                  }
                  else {
                    mBeta[j]     = 1.0;
                    mBetaStar[j] = 0.0;
                    mPUTag[j]    = 1.0;
                  }
                  mPhi[j] = mEvent->pfjet(j).phi();
                  mEta[j] = mEvent->pfjet(j).eta();
                  mMass[j]= mEvent->pfjet(j).mass() * mEvent->pfjet(j).cor();
                  mCHF[j] = mEvent->pfjet(j).chf();
                  mNHF[j] = mEvent->pfjet(j).nhf();
                  mPHF[j] = mEvent->pfjet(j).phf();
                  mELF[j] = mEvent->pfjet(j).elf();
                  mMUF[j] = mEvent->pfjet(j).muf();
                }// pfjet loop
                mOutTree->Fill();
                counter_id++;
              }// id cuts
            }// at least 8 pfjets
          }// if pv
        }// if hlt
      }// trigger loop
    }// tree loop
  }// file loop
  if (Ntrig > 0) {
    cout<<"Events passing the trigger:       "<<counter_hlt<<endl;
  } 
  cout<<"Events passing the PV:            "<<counter_pv<<endl;
  cout<<"Events passing the HCAL NOISE:    "<<counter_hcal<<endl;  
  cout<<"Events passing id & Kin cuts:     "<<counter_id<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
vector<int> MultijetSearchTree::findIndices(const vector<int>& v, int rank)
{
  vector<int> indices;
  for(int i=0;i<2*rank;i++) {
    bool found(false);
    for(int j=0;j<rank;j++) {
      if (i == v[j]) {
        found = true;
        continue;
      }
    }
    if (!found) {    
      indices.push_back(i);
    }
  }
  //cout<<indices.size()<<" "<<indices[0]<<indices[1]<<indices[2]<<indices[3]<<endl;
  return indices;
}
//////////////////////////////////////////////////////////////////////////////////////////
MultijetSearchTree::~MultijetSearchTree() 
{
}

DEFINE_FWK_MODULE(MultijetSearchTree);
