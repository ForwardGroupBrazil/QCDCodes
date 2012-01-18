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

#include "KKousour/MultiJetAnalysis/plugins/PatMultijetSearchTree.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

using namespace std;
using namespace reco;

PatMultijetSearchTree::PatMultijetSearchTree(edm::ParameterSet const& cfg) 
{
  srcJets_     = cfg.getParameter<edm::InputTag>             ("jets");
  srcMET_      = cfg.getParameter<edm::InputTag>             ("met");
  srcBeta_     = cfg.getParameter<string>                    ("beta"); 
  mEtaMax      = cfg.getParameter<double>                    ("etaMAX");
  mPtMin       = cfg.getParameter<double>                    ("ptMIN");
  mBetaMax     = cfg.getParameter<double>                    ("betaMAX");
  vector<double> vLumi(0),vBnd(0);
  mIsMC        = cfg.getUntrackedParameter<bool>             ("isMC",false);
  mPtHatLumi   = cfg.getUntrackedParameter<vector<double> >  ("ptHatLumi",vLumi);
  mPtHatBnd    = cfg.getUntrackedParameter<vector<double> >  ("ptHatBnd",vBnd);
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetSearchTree::beginJob() 
{
  mOutTree = fs->make<TTree>("events","events");
  mOutTree->Branch("runNo"       ,&mRun         ,"mRun/I");
  mOutTree->Branch("evtNo"       ,&mEvt         ,"mEvt/I");
  mOutTree->Branch("npv"         ,&mNPV         ,"mNPV/I");
  mOutTree->Branch("m4jIndex"    ,&mM4JIndex    ,"mM4JIndex[2]/I");
  mOutTree->Branch("m2jIndex"    ,&mM2JIndex    ,"mM2JIndex[2][2]/I");
  mOutTree->Branch("metSig"      ,&mMetSig      ,"mMetSig/F");
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
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetSearchTree::endJob() 
{

}
//////////////////////////////////////////////////////////////////////////////////////////
int PatMultijetSearchTree::getBin(double x, const std::vector<double>& boundaries)
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
void PatMultijetSearchTree::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  //----- MC -------------------
  double wt(1.0);
  /*
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
  */
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(srcJets_,jets);
  edm::View<pat::Jet> pat_jets = *jets;

  edm::Handle<edm::View<MET> >  met;
  iEvent.getByLabel(srcMET_,met);

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("goodOfflinePrimaryVertices",recVtxs);

  bool cut_vtx = (recVtxs->size() > 0);
  bool cut_njets = (pat_jets.size() > 7);

  if (cut_vtx && cut_njets) {
    double HT(0),HT4J_1(0),HT4J_2(0),M4J(0),M4J_1(0),M4J_2(0),cosThetaStar(0);
    bool cutID(true),cut_eta(true),cut_pt(true),cut_pu(true);
    LorentzVector P4[8],P48J(0,0,0,0);
    int N(0);
    for(edm::View<pat::Jet>::const_iterator ijet = pat_jets.begin();ijet != pat_jets.end() && N < 8; ++ijet) { 
      bool id(false);
      double chf = ijet->chargedHadronEnergyFraction();
      double nhf = ijet->neutralHadronEnergyFraction() + ijet->HFHadronEnergyFraction();
      double phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double muf = ijet->muonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      int chm = ijet->chargedHadronMultiplicity();
      int npr = ijet->chargedMultiplicity() + ijet->neutralMultiplicity();
      id = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(ijet->eta())<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || fabs(ijet->eta())>2.4));
      double beta = ijet->userFloat(srcBeta_);
      cut_pu *= (beta < mBetaMax);
      cutID *= id;
      cut_pt  *= (ijet->pt() > mPtMin);
      cut_eta *= (fabs(ijet->eta()) < mEtaMax);
      HT += ijet->pt();
      P4[N] = ijet->p4();
      P48J += P4[N];
      mPt[N]  = ijet->pt();
      mPhi[N] = ijet->phi();
      mEta[N] = ijet->eta();
      mMass[N]= ijet->mass();
      mBeta[N]= beta;
      mCHF[N] = chf;
      mNHF[N] = nhf;
      mPHF[N] = phf;
      mELF[N] = elf;
      mMUF[N] = muf;
      N++;
    }// jet loop
    if (cutID && cut_pu && cut_eta && cut_pt) {
      mM8J    = P48J.mass();
      mHT     = HT;
      mMetSig = (*met)[0].et()/(*met)[0].sumEt();
      mRun    = iEvent.id().run();
      mEvt    = iEvent.id().event();
      mNPV    = recVtxs->size();
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
      mM4JBalance = 2*fabs(mM4J[0]-mM4J[1])/(mM4J[0]+mM4J[1]);
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
      mOutTree->Fill();
    }// if jet cuts
  }// if vtx and jet multi
}
//////////////////////////////////////////////////////////////////////////////////////////
vector<int> PatMultijetSearchTree::findIndices(const vector<int>& v, int rank)
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
PatMultijetSearchTree::~PatMultijetSearchTree() 
{
}

DEFINE_FWK_MODULE(PatMultijetSearchTree);
