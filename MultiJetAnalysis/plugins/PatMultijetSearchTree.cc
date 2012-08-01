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
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

using namespace std;
using namespace reco;

PatMultijetSearchTree::PatMultijetSearchTree(edm::ParameterSet const& cfg) 
{
  srcJets_    = cfg.getParameter<edm::InputTag>          ("jets");
  srcMET_     = cfg.getParameter<edm::InputTag>          ("met");
  srcRho_     = cfg.getParameter<edm::InputTag>          ("rho");
  srcGenJets_ = cfg.getUntrackedParameter<edm::InputTag> ("genjets",edm::InputTag("ak5GenJets"));
  srcPU_      = cfg.getUntrackedParameter<string>        ("pu","");
  isMC_       = cfg.getUntrackedParameter<bool>          ("isMC",false);
  etaMax_     = cfg.getParameter<double>                 ("etaMAX");
  ptMin_      = cfg.getParameter<double>                 ("ptMIN");
  betaMin_    = cfg.getParameter<double>                 ("betaMIN");
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetSearchTree::beginJob() 
{
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"       ,&run_         ,"run_/I");
  outTree_->Branch("evtNo"       ,&evt_         ,"evt_/I");
  outTree_->Branch("nvtx"        ,&nVtx_        ,"nVtx_/I");
  outTree_->Branch("index2J"     ,&index2J_     ,"index2J_[4][2]/I");
  outTree_->Branch("index4J"     ,&index4J_     ,"index4J_[2][4]/I");
  outTree_->Branch("rho"         ,&rho_         ,"rho_/F");
  outTree_->Branch("metSig"      ,&metSig_      ,"metSig_/F");
  outTree_->Branch("dphi4j"      ,&dPhi4j_      ,"dPhi4j_/F");
  outTree_->Branch("ht"          ,&ht_          ,"ht_/F");
  outTree_->Branch("m8j"         ,&m8j_         ,"m8j_/F");
  outTree_->Branch("m4jAve"      ,&m4jAve_      ,"m4jAve_/F");
  outTree_->Branch("m2jAve"      ,&m2jAve_      ,"m2jAve_/F");
  outTree_->Branch("m2jEst"      ,&m2jEst_      ,"m2jEst_/F");
  outTree_->Branch("m2jSig"      ,&m2jSigma_    ,"m2jSigma_/F");
  outTree_->Branch("cosThetaStar",&cosThetaStar_,"cosThetaStar_/F");
  outTree_->Branch("m4jBalance"  ,&m4jBalance_  ,"m4jBalance_/F");
  outTree_->Branch("m2j"         ,&m2j_         ,"m2j_[4]/F");
  outTree_->Branch("m4j"         ,&m4j_         ,"m4j_[2]_/F");
  outTree_->Branch("ht4j"        ,&ht4j_        ,"ht4j_[2]/F");
  outTree_->Branch("eta4j"       ,&eta4j_       ,"eta4j_[2]/F");
  outTree_->Branch("pt4j"        ,&pt4j_        ,"pt4j_[2]/F");
  outTree_->Branch("dR2jAll"     ,&dR2jAll_     ,"dR2jAll_[28]/F");
  pt_   = new std::vector<float>;
  eta_  = new std::vector<float>;
  phi_  = new std::vector<float>;
  mass_ = new std::vector<float>;
  jec_  = new std::vector<float>;
  unc_  = new std::vector<float>;
  beta_ = new std::vector<float>;
  chf_  = new std::vector<float>;
  nhf_  = new std::vector<float>;
  phf_  = new std::vector<float>;
  elf_  = new std::vector<float>;
  muf_  = new std::vector<float>;
  outTree_->Branch("pt"          ,"vector<float>" ,&pt_);
  outTree_->Branch("jec"         ,"vector<float>" ,&jec_);
  outTree_->Branch("unc"         ,"vector<float>" ,&unc_);
  outTree_->Branch("beta"        ,"vector<float>" ,&beta_);
  outTree_->Branch("eta"         ,"vector<float>" ,&eta_);
  outTree_->Branch("phi"         ,"vector<float>" ,&phi_);
  outTree_->Branch("mass"        ,"vector<float>" ,&mass_);
  outTree_->Branch("chf"         ,"vector<float>" ,&chf_);
  outTree_->Branch("nhf"         ,"vector<float>" ,&nhf_);
  outTree_->Branch("phf"         ,"vector<float>" ,&phf_);
  outTree_->Branch("muf"         ,"vector<float>" ,&muf_);
  outTree_->Branch("elf"         ,"vector<float>" ,&elf_);
  //------------------- MC ---------------------------------
  if (isMC_) {
    outTree_->Branch("pu"             ,&simPU_         ,"simPU_/I");
    outTree_->Branch("m4jAveGEN"      ,&m4jAveGEN_     ,"m4jAveGEN_/F");
    outTree_->Branch("m2jAveGEN"      ,&m2jAveGEN_     ,"m2jAveGEN_/F");
    outTree_->Branch("m2jEstGEN"      ,&m2jEstGEN_     ,"m2jEstGEN_/F");
    outTree_->Branch("m4jAveParton"   ,&m4jAveParton_  ,"m4jAveParton_/F");
    outTree_->Branch("m2jAveParton"   ,&m2jAveParton_  ,"m2jAveParton_/F");
    outTree_->Branch("m2jEstParton"   ,&m2jEstParton_  ,"m2jEstParton_/F");
    partonId_  = new std::vector<int>;
    partonSt_  = new std::vector<int>;
    partonPt_  = new std::vector<float>;
    partonEta_ = new std::vector<float>;
    partonPhi_ = new std::vector<float>;
    partonE_   = new std::vector<float>;
    outTree_->Branch("partonId"       ,"vector<int>"   ,&partonId_);
    outTree_->Branch("partonSt"       ,"vector<int>"   ,&partonSt_);
    outTree_->Branch("partonPt"       ,"vector<float>" ,&partonPt_);
    outTree_->Branch("partonEta"      ,"vector<float>" ,&partonEta_);
    outTree_->Branch("partonPhi"      ,"vector<float>" ,&partonPhi_);
    outTree_->Branch("partonE"        ,"vector<float>" ,&partonE_);
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetSearchTree::endJob() 
{
  if (isMC_) {
    delete partonSt_;
    delete partonId_;
    delete partonPt_;
    delete partonEta_;
    delete partonPhi_;
    delete partonE_;
  }
  delete pt_;
  delete eta_;
  delete phi_;
  delete mass_;
  delete jec_;
  delete unc_;
  delete chf_;
  delete nhf_;
  delete phf_;
  delete elf_;
  delete muf_;
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetSearchTree::initialize()
{
  run_          = -999;
  evt_          = -999;
  nVtx_         = -999;
  rho_          = -999;
  metSig_       = -999;
  dPhi4j_       = -999;
  ht_           = -999;
  m8j_          = -999;
  m4jAve_       = -999;
  m2jAve_       = -999;
  m2jEst_       = -999;
  m2jSigma_     = -999;
  cosThetaStar_ = -999;
  m4jBalance_   = -999;
  for(int i=0;i<2;i++) {
    m4j_[i]   = -999;
    ht4j_[i]  = -999;
    eta4j_[i] = -999;
    pt4j_[i]  = -999;
    for(int j=0;j<4;j++) {
      index2J_[j][i] = -999;
      index4J_[i][j] = -999; 
    }
  }
  for(int i=0;i<4;i++) {
    m2j_[i] = -999;
  }
  for(int i=0;i<28;i++) {
    dR2jAll_[i] = -999;
  }
  pt_->clear();
  eta_->clear();
  phi_->clear();
  mass_->clear();
  chf_->clear();
  nhf_->clear();
  phf_->clear();
  elf_->clear();
  muf_->clear();
  jec_->clear();
  unc_->clear();
  //----- MC -------
  if (isMC_) {
    simPU_        = -999;
    m2jAveGEN_    = -999;
    m2jEstGEN_    = -999;
    m4jAveGEN_    = -999;
    m2jAveParton_ = -999;
    m2jEstParton_ = -999;
    m4jAveParton_ = -999;
    partonSt_ ->clear();
    partonId_ ->clear();
    partonPt_ ->clear();
    partonEta_->clear();
    partonPhi_->clear();
    partonE_  ->clear();
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetSearchTree::getDoublets(const vector<LorentzVector>& P4, float& m2jAve, float& m2jEst, float& m2jSigma, float m2j[], int index2J[][2],int index[])
{
  //----- construct all possible 4 doublet configurations
  unsigned int NJETS = P4.size();
  double sigMIN(10000),dmMIN(10000);
  for(unsigned k5=5;k5<NJETS-2;k5++) {
    for(unsigned k6=k5+1;k6<NJETS-1;k6++) {
      for(unsigned k7=k6+1;k7<NJETS;k7++) {
        vector<LorentzVector> tmpP4; 
        int trueIndex[8];
        for(int j=0;j<5;j++) {
          tmpP4.push_back(P4[j]);
          trueIndex[j] = j;
        }
        tmpP4.push_back(P4[k5]);
        tmpP4.push_back(P4[k6]);
        tmpP4.push_back(P4[k7]);
        trueIndex[5] = k5;
        trueIndex[6] = k6;
        trueIndex[7] = k7;
        LorentzVector pairP4[4][2];
        int pairInd[4][2];
        for(int j1=0;j1<8;j1++) {
          for(int j2=j1+1;j2<8;j2++) {
            if (j2==j1) continue;
            for(int j3=0;j3<8;j3++) {
              if (j3==j2 || j3==j1) continue;
              for(int j4=j3+1;j4<8;j4++) {
                if (j4==j3 || j4==j2 || j4==j1) continue;
                for(int j5=0;j5<8;j5++) { 
                  if (j5==j4 || j5==j3 || j5==j2 || j5==j1) continue;
                  for(int j6=j5+1;j6<8;j6++) {
                    if (j6==j5 || j6==j4 || j6==j3 || j6==j2 || j6==j1) continue;
                    for(int j7=0;j7<8;j7++) {
                      if (j7==j6 || j7==j5 || j7==j4 || j7==j3 || j7==j2 || j7==j1) continue;
                      for(int j8=j7+1;j8<8;j8++) {
                        if (j8==j7 || j8==j6 || j8==j5 || j8==j4 || j8==j3 || j8==j2 || j8==j1) continue;
                        pairP4[0][0] = tmpP4[j1];
                        pairP4[0][1] = tmpP4[j2];
                        pairP4[1][0] = tmpP4[j3];
                        pairP4[1][1] = tmpP4[j4];
                        pairP4[2][0] = tmpP4[j5];
                        pairP4[2][1] = tmpP4[j6];
                        pairP4[3][0] = tmpP4[j7];
                        pairP4[3][1] = tmpP4[j8];
                        pairInd[0][0] = trueIndex[j1];
                        pairInd[0][1] = trueIndex[j2];
                        pairInd[1][0] = trueIndex[j3];
                        pairInd[1][1] = trueIndex[j4];
                        pairInd[2][0] = trueIndex[j5];
                        pairInd[2][1] = trueIndex[j6];
                        pairInd[3][0] = trueIndex[j7];
                        pairInd[3][1] = trueIndex[j8]; 
                        double sum(0.0),sum2(0.0),M2J[4];
                        for(int ip=0;ip<4;ip++) {
                          M2J[ip] = (pairP4[ip][0]+pairP4[ip][1]).mass();
                          sum  += M2J[ip];
                          sum2 += M2J[ip]*M2J[ip];
                        }
                        double ave = 0.25*sum;
                        double sig = sqrt((sum2-4*ave*ave)/3.0); 
                        double min_pair_diff(10000);
                        double tmp_m2jEst(0.0);
                        for(int ip1=0;ip1<4;ip1++) {
                          for(int ip2=ip1+1;ip2<4;ip2++) {
                            double dm = fabs(M2J[ip1]-M2J[ip2]);
                            if (dm < min_pair_diff) {
                              min_pair_diff = dm;
                              tmp_m2jEst = 0.5*(M2J[ip1]+M2J[ip2]);
                            }
                          }
                        }
                        if (min_pair_diff < dmMIN) {
                          dmMIN = min_pair_diff;
                          m2jEst = tmp_m2jEst;
                          for(int i=0;i<8;i++) {
                            index[i] = trueIndex[i];
                          }
                          for(int ip=0;ip<4;ip++) {
                            m2j[ip] = M2J[ip];
                            index2J[ip][0] = pairInd[ip][0];
                            index2J[ip][1] = pairInd[ip][1];
                          } 
                        } 
                        //cout<<j1<<" "<<j2<<" "<<j3<<" "<<j4<<" "<<j5<<" "<<j6<<" "<<j7<<" "<<j8<<" "<<ave<<" "<<sig<<endl;
                        if (sig/ave < sigMIN) {
                          m2jAve   = ave;
                          m2jSigma = sig/ave;
                          /*
                          for(int i=0;i<8;i++) {
                            index[i] = trueIndex[i];
                          }
                          for(int ip=0;ip<4;ip++) {
                            m2j[ip] = M2J[ip];
                            index2J[ip][0] = pairInd[ip][0];
                            index2J[ip][1] = pairInd[ip][1];
                          }
                          */
                          sigMIN = sig/ave;
                        }       
                      }// j8
                    }// j7
                  }// j6
                }// j5 
              }//j4 
            }// j3
          }// j2  
        }//j1
      }// k7
    }// k6
  }// k5
  /*
  cout<<"Found best combination from "<<NJETS<<" jets: ";
  for(int j=0;j<8;j++)
    cout<<index[j]<<" ";
  cout<<endl;
  */
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetSearchTree::getQuartets(const std::vector<LorentzVector>& P4, LorentzVector quartetP4[], float& m4jAve, float& m4jBalance, float m4j[], int index4J[][4], int index2J[][2])
{
  //---- loop over the optimum set of pairs --------
  vector<int> v[2];
  double d,dmin(10000),quartetM[2];
  int q1[2],q2[2],q3[2],q4[2];
  for(int jj1=0;jj1<1;jj1++){
    for(int jj2=jj1+1;jj2<4;jj2++){
      v[0].clear();
      v[1].clear();
      v[0].push_back(jj1);
      v[0].push_back(jj2);
      // find the complementary indices
      v[1] = findIndices(v[0],2);
      for(int iq=0;iq<2;iq++) {
        q1[iq] = index2J[v[iq][0]][0];
        q2[iq] = index2J[v[iq][0]][1];
        q3[iq] = index2J[v[iq][1]][0];  
        q4[iq] = index2J[v[iq][1]][1];
        quartetP4[iq] = P4[q1[iq]]+P4[q2[iq]]+P4[q3[iq]]+P4[q4[iq]];
        quartetM[iq] = quartetP4[iq].mass(); 
      } 
      d = abs(quartetM[0]-quartetM[1]);
      if (d < dmin) {
        m4jAve       = 0.5*(quartetM[0]+quartetM[1]);
        m4jBalance   = fabs(quartetM[0]-quartetM[1])/m4jAve;
        for(int iq=0;iq<2;iq++) {
          m4j[iq]   = quartetP4[iq].mass();
          index4J[iq][0] = q1[iq];
          index4J[iq][1] = q2[iq];
          index4J[iq][2] = q3[iq];
          index4J[iq][3] = q4[iq];
        }
        dmin = d;
      } 
    }
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

  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("goodOfflinePrimaryVertices",recVtxs);

  edm::Handle<GenEventInfoProduct> hEventInfo;
  edm::Handle<GenJetCollection> genjets;
  edm::Handle<GenParticleCollection> genParticles;
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;

  initialize();

  int intpu(0);
  if (!iEvent.isRealData()) {
    if (srcPU_ != "") {
      iEvent.getByLabel(srcPU_,PupInfo); 
      for(vector<PileupSummaryInfo>::const_iterator ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu) {
        if (ipu->getBunchCrossing() == 0) {
          intpu += ipu->getPU_NumInteractions();
        } 
      }
    }
    //---------- genjets ------------------
    iEvent.getByLabel(srcGenJets_,genjets);
    vector<LorentzVector> genP4;
    for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
      if (i_gen->pt() > ptMin_ && fabs(i_gen->eta()) < etaMax_)
        genP4.push_back(i_gen->p4());
    }
    if (genP4.size() > 7) {  
      float m2jGEN[4],m2jSigmaGEN,m4jBalanceGEN,m4jGEN[4];
      int index2JGEN[4][2],index4JGEN[2][4],indexGEN[8];
      getDoublets(genP4,m2jAveGEN_,m2jEstGEN_,m2jSigmaGEN,m2jGEN,index2JGEN,indexGEN);
      LorentzVector quartetP4GEN[2];
      getQuartets(genP4,quartetP4GEN,m4jAveGEN_,m4jBalanceGEN,m4jGEN,index4JGEN,index2JGEN);
    }
    //---------- partons ------------------
    iEvent.getByLabel("genParticles", genParticles);
    vector<LorentzVector> partonP4;
    for(unsigned ip = 0; ip < genParticles->size(); ++ ip) {
      const GenParticle &p = (*genParticles)[ip];
      if (p.status() != 3) continue;
      int ndst3(0);
      for(unsigned k = 0; k < p.numberOfDaughters(); k++) {
        if (p.daughter(k)->status() == 3) {
          ndst3++; 
        }
      }
      if (ndst3 > 0) continue;
      partonId_ ->push_back(p.pdgId());
      partonSt_ ->push_back(p.status()); 
      partonPt_ ->push_back(p.pt());
      partonEta_->push_back(p.eta());
      partonPhi_->push_back(p.phi());
      partonE_  ->push_back(p.energy()); 
      partonP4.push_back(p.p4());
    }// parton loop
    if (partonP4.size() > 7) {
      float m2jParton[4],m2jSigmaParton,m4jBalanceParton,m4jParton[2];
      int index2JParton[4][2],index4JParton[2][4],indexParton[8];
      getDoublets(partonP4,m2jAveParton_,m2jEstParton_,m2jSigmaParton,m2jParton,index2JParton,indexParton);
      LorentzVector quartetP4Parton[2];
      getQuartets(partonP4,quartetP4Parton,m4jAveParton_,m4jBalanceParton,m4jParton,index4JParton,index2JParton);
    }
  } 
  bool cut_vtx = (recVtxs->size() > 0);
  bool cut_njets = (pat_jets.size() > 7);

  if (cut_vtx && cut_njets) {
    bool cutID(true),cut_eta(true),cut_pt(true),cut_pu(true);
    vector<LorentzVector> P4(0); 
    LorentzVector P48J(0,0,0,0);
    int N(0);
    double HT(0.0);
    for(edm::View<pat::Jet>::const_iterator ijet = pat_jets.begin();ijet != pat_jets.end(); ++ijet) { 
      bool id(false);
      double chf = ijet->chargedHadronEnergyFraction();
      double nhf = ijet->neutralHadronEnergyFraction() + ijet->HFHadronEnergyFraction();
      double phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double muf = ijet->muonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      int chm = ijet->chargedHadronMultiplicity();
      int npr = ijet->chargedMultiplicity() + ijet->neutralMultiplicity();
      id = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(ijet->eta())<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || fabs(ijet->eta())>2.4));
      double beta = ijet->userFloat("beta");
      cut_pu  = (beta > betaMin_);
      cutID   = id;
      cut_pt  = (ijet->pt() > ptMin_);
      cut_eta = (fabs(ijet->eta()) < etaMax_);
      if (cut_pu && cutID && cut_pt && cut_eta) {
        P4.push_back(ijet->p4());
        pt_  ->push_back(ijet->pt());
        phi_ ->push_back(ijet->phi());
        eta_ ->push_back(ijet->eta());
        mass_->push_back(ijet->mass());
        beta_->push_back(beta);
        chf_ ->push_back(chf);
        nhf_ ->push_back(nhf);
        phf_ ->push_back(phf);
        elf_ ->push_back(elf);
        muf_ ->push_back(muf);
        jec_ ->push_back(1./ijet->jecFactor(0));
        unc_ ->push_back(ijet->userFloat("jecUnc"));
        N++;
      }
    }// jet loop
    if (N > 7) {
      rho_    = *rho;
      metSig_ = (*met)[0].et()/(*met)[0].sumEt();
      run_    = iEvent.id().run();
      evt_    = iEvent.id().event();
      nVtx_   = recVtxs->size();
      simPU_  = intpu;
      getDoublets(P4,m2jAve_,m2jEst_,m2jSigma_,m2j_,index2J_,index_);
      for(int j=0;j<8;j++) {
        HT += P4[index_[j]].pt();  
        P48J += P4[index_[j]];
      }
      ht_  = HT;
      m8j_ = P48J.mass();
      LorentzVector quartetP4[2];
      getQuartets(P4,quartetP4,m4jAve_,m4jBalance_,m4j_,index4J_,index2J_);     
      cosThetaStar_ = tanh(0.5*(quartetP4[0].Rapidity()-quartetP4[1].Rapidity()));
      dPhi4j_       = abs(deltaPhi(quartetP4[0].phi(),quartetP4[1].phi()));
      for(int iq=0;iq<2;iq++) {
        ht4j_[iq]  = P4[index4J_[iq][0]].pt()+P4[index4J_[iq][1]].pt()+P4[index4J_[iq][2]].pt()+P4[index4J_[iq][3]].pt();
        eta4j_[iq] = quartetP4[iq].eta();
        pt4j_[iq]  = quartetP4[iq].pt();
      }  
      //----- 2J combinatorics -------- 
      int jj(0);
      for(int j1=0;j1<8;j1++) {
        for(int j2=j1+1;j2<8;j2++) {
          dR2jAll_[jj] = deltaR(P4[j1],P4[j2]);
          jj++;
        }
      } 
      outTree_->Fill();
    }// if at least 8 good jets
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
