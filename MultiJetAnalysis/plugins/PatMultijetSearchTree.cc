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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

using namespace std;
using namespace reco;

PatMultijetSearchTree::PatMultijetSearchTree(edm::ParameterSet const& cfg) 
{
  srcJets_ = cfg.getParameter<edm::InputTag> ("jets");
  srcMET_  = cfg.getParameter<edm::InputTag> ("met");
  srcBeta_ = cfg.getParameter<string>        ("beta"); 
  etaMax_  = cfg.getParameter<double>        ("etaMAX");
  ptMin_   = cfg.getParameter<double>        ("ptMIN");
  betaMax_ = cfg.getParameter<double>        ("betaMAX");
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetSearchTree::beginJob() 
{
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"       ,&run_         ,"run_/I");
  outTree_->Branch("evtNo"       ,&evt_         ,"evt_/I");
  outTree_->Branch("nvtx"        ,&nVtx_        ,"nVtx_/I");
  outTree_->Branch("pu"          ,&simPU_       ,"simPU_/I");
  outTree_->Branch("m4jIndex"    ,&m4JIndex_    ,"m4JIndex_[2]/I");
  outTree_->Branch("m2jIndex"    ,&m2JIndex_    ,"m2JIndex_[2][2]/I");
  outTree_->Branch("metSig"      ,&metSig_      ,"metSig_/F");
  outTree_->Branch("dphi4j"      ,&dPhi4J_      ,"dPhi4J_/F");
  outTree_->Branch("ht"          ,&ht_          ,"ht_/F");
  outTree_->Branch("m8j"         ,&m8J_         ,"m8J_/F");
  outTree_->Branch("cosThetaStar",&cosThetaStar_,"cosThetaStar_/F");
  outTree_->Branch("m4jBalance"  ,&m4JBalance_  ,"m4JBalance_/F");
  outTree_->Branch("m2jBalance"  ,&m2JBalance_  ,"m2JBalance_[2]/F"); 
  outTree_->Branch("m4j"         ,&m4J_         ,"m4J_[2]_/F");
  outTree_->Branch("ht4j"        ,&ht4J_        ,"ht4J_[2]/F");
  outTree_->Branch("eta4j"       ,&eta4J_       ,"eta4J_[2]/F");
  outTree_->Branch("pt4j"        ,&pt4J_        ,"pt4J_[2]/F");
  outTree_->Branch("m2j"         ,&m2J_         ,"m2J_[2][2]/F");
  outTree_->Branch("dR2jAll"     ,&dR2JAll_     ,"dR2JAll_[28]/F");
  outTree_->Branch("pt"          ,&pt_          ,"pt_[8]/F");
  outTree_->Branch("beta"        ,&beta_        ,"beta_[8]/F");
  outTree_->Branch("eta"         ,&eta_         ,"eta_[8]/F");
  outTree_->Branch("phi"         ,&phi_         ,"phi_[8]/F");
  outTree_->Branch("mass"        ,&mass_        ,"mass_[8]/F");
  outTree_->Branch("chf"         ,&chf_         ,"chf_[8]/F");
  outTree_->Branch("nhf"         ,&nhf_         ,"nhf_[8]/F");
  outTree_->Branch("phf"         ,&phf_         ,"phf_[8]/F");
  outTree_->Branch("muf"         ,&muf_         ,"muf_[8]/F");
  outTree_->Branch("elf"         ,&elf_         ,"elf_[8]/F");
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetSearchTree::endJob() 
{

}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetSearchTree::initialize()
{
  run_          = -999;
  evt_          = -999;
  nVtx_         = -999;
  simPU_        = -999;
  metSig_       = -999;
  dPhi4J_       = -999;
  ht_           = -999;
  m8J_          = -999;
  cosThetaStar_ = -999;
  m4JBalance_   = -999;
  for(int i=0;i<2;i++) {
    m4JIndex_[i]   = -999;
    m2JBalance_[i] = -999;
    m4J_[i]        = -999;
    ht4J_[i]       = -999;
    eta4J_[i]      = -999;
    pt4J_[i]       = -999;
    for(int j=0;j<2;j++) {
      m2JIndex_[i][j] = -999;
      m2J_[i][j]      = -999;
    }
  }
  for(int i=0;i<28;i++) {
    dR2JAll_[i] = -999;
  }
  for(int i=0;i<8;i++) {
    pt_[i]   = -999;
    eta_[i]  = -999;
    phi_[i]  = -999;
    mass_[i] = -999;
    chf_[i]  = -999;
    nhf_[i]  = -999;
    phf_[i]  = -999;
    elf_[i]  = -999;
    muf_[i]  = -999;
  }  
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetSearchTree::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(srcJets_,jets);
  edm::View<pat::Jet> pat_jets = *jets;

  edm::Handle<edm::View<MET> >  met;
  iEvent.getByLabel(srcMET_,met);

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("goodOfflinePrimaryVertices",recVtxs);

  edm::Handle<GenEventInfoProduct> hEventInfo;
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;

  int intpu(0);
  if (!iEvent.isRealData()) {
    iEvent.getByLabel("addPileupInfo",PupInfo);
    for(vector<PileupSummaryInfo>::const_iterator ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu) {
      if (ipu->getBunchCrossing() == 0)
        intpu += ipu->getPU_NumInteractions(); 
    }
  } 
  bool cut_vtx = (recVtxs->size() > 0);
  bool cut_njets = (pat_jets.size() > 7);

  if (cut_vtx && cut_njets) {
    double HT(0),HT4J_1(0),HT4J_2(0),M4J_1(0),M4J_2(0),cosThetaStar(0);
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
      cut_pu *= (beta < betaMax_);
      cutID *= id;
      cut_pt  *= (ijet->pt() > ptMin_);
      cut_eta *= (fabs(ijet->eta()) < etaMax_);
      HT += ijet->pt();
      P4[N] = ijet->p4();
      P48J += P4[N];
      pt_[N]  = ijet->pt();
      phi_[N] = ijet->phi();
      eta_[N] = ijet->eta();
      mass_[N]= ijet->mass();
      beta_[N]= beta;
      chf_[N] = chf;
      nhf_[N] = nhf;
      phf_[N] = phf;
      elf_[N] = elf;
      muf_[N] = muf;
      N++;
    }// jet loop
    if (cutID && cut_pu && cut_eta && cut_pt) {
      m8J_    = P48J.mass();
      ht_     = HT;
      metSig_ = (*met)[0].et()/(*met)[0].sumEt();
      run_    = iEvent.id().run();
      evt_    = iEvent.id().event();
      nVtx_   = recVtxs->size();
      simPU_  = intpu;
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
                m4J_[0] = M4J_1;
                m4J_[1] = M4J_2;
                eta4J_[0] = P44J_1.eta();
                eta4J_[1] = P44J_2.eta();
                pt4J_[0] = P44J_1.pt();
                pt4J_[1] = P44J_2.pt();
                ht4J_[0] = HT4J_1;
                ht4J_[1] = HT4J_2;
                cosThetaStar_ = cosThetaStar;
                dPhi4J_ = deltaPhi(P44J_1.phi(),P44J_2.phi());  
                m4JIndex_[0] = (v1[0]+1)*1000+(v1[1]+1)*100+(v1[2]+1)*10+(v1[3]+1);
                m4JIndex_[1] = (v2[0]+1)*1000+(v2[1]+1)*100+(v2[2]+1)*10+(v2[3]+1);
                v4J[0] = v1;
                v4J[1] = v2;
                dmin = d;
              } 
            }
          }
        }
      }
      m4JBalance_ = 2*fabs(m4J_[0]-m4J_[1])/(m4J_[0]+m4J_[1]);
      //----- 2J combinatorics -------- 
      int jj(0);
      for(int j1=0;j1<8;j1++) {
        for(int j2=j1+1;j2<8;j2++) {
          dR2JAll_[jj] = deltaR(P4[j1],P4[j2]);
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
              m2J_[k][0] = M2J[0];
              m2J_[k][1] = M2J[1];
              m2JBalance_[k] = 2*fabs(m2J_[k][0]-m2J_[k][1])/(m2J_[k][0]+m2J_[k][1]);
              m2JIndex_[k][0] = (v4J[k][v1[0]]+1)*10+(v4J[k][v1[1]]+1);
              m2JIndex_[k][1] = (v4J[k][v2[0]]+1)*10+(v4J[k][v2[1]]+1);
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
      outTree_->Fill();
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
  return indices;
}
//////////////////////////////////////////////////////////////////////////////////////////
PatMultijetSearchTree::~PatMultijetSearchTree() 
{
}

DEFINE_FWK_MODULE(PatMultijetSearchTree);
