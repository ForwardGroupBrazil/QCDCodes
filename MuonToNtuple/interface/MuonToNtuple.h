#ifndef MuonToNtuple_H
#define MuonToNtuple_H

/*
1. muon
- 3 vector momentum
- tag and probe info (for example, the muon passes the requirement of
"tag", tag = 1. If not, tag = 0. Same to the probe)

2. Jet
- kinematic variables (pt, 3 v momentum, eta, phi...) for each jet

3. global event information
- luminosity, run number, event number

4. Primary vertex (PV) information
- number of primary vertices
- PVz (z coordinate), # of tracks for the PV

5. track-matching information for the muon
- track z, phi, eta, pt,
- # of hits in trackers for each track

6. trigger efficiency
- when it's available

7. electron track information
- to know the electron track is sharing with muon track-matching (it can
be a kind of background)

8. for cosmic ray bg
- timing information ...

9. various topological variables (may be useful for the physics
background suppression such as ttbar, Z->tautau)
- H_T, Spericity, Centrality, Aplanarity ...
- I have a source code for these variables and I think we can implement
it to your code easily (or as another class).
*/
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "DataFormats/Common/interface/RefToBaseVector.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h"


#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TMath.h>

//
// class decleration
//

using namespace std;
namespace reco { class CandCommonVertexFitterBase; class VertexCompositeCandidate; class CandCommonVertexFitter; }

class MuonToNtuple : public edm::EDAnalyzer {
   public:
      explicit MuonToNtuple(const edm::ParameterSet&);
      ~MuonToNtuple();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;


      edm::RefToBaseVector<pat::Muon> firstMuon(const edm::View<pat::Muon>&) const;

      edm::RefToBaseVector<pat::Muon> secondMuon(const edm::View<pat::Muon>&, const pat::Muon*) const;

      float angleBetween(const pat::Muon& lhs, const pat::Muon& rhs) const;
      float angleBetween(const reco::Track& lhs, const reco::Track& rhs) const;
      float angleBetween(const reco::GenParticle& lhs, const reco::GenParticle& rhs) const;

      float deltaEtaPhi(const pat::Muon& lhs, const pat::Muon& rhs) const;

      int motherId(const reco::GenParticle& par) const;
      int motherId(const reco::GenParticleRef par) const;

      double TrackSumPtrInCone03( const pat::Muon& muon, edm::Handle<edm::View<reco::Track> > tracks );

      std::string theRootFileName;
      edm::InputTag theMuonLabel;
      edm::InputTag theMETLabel;
      edm::InputTag theJetLabel;
      double theCrossSection;
      double theFilterEfficiency;
      double theTotalNevents;
      double theBDiscriminant;
      bool isMC;

      TFile *theFile;
      TTree *DiMuonTree;


      int nEvt;

      //FILL a tree
      static const int MPSIZE = 1000;

      // Invariant Mass distribution of SS(OS) di-muon events GG (2 global)
      // GlbTree
      int runNum;
      int evtNum;
      int nPair;
      int GENnPair;
      double weight;
      double MET;
      int Njets;
      int Nmuons;
      int Nbtagged;
      int NbtaggedCloseMuon;

      // PV
      int PVtrackSize;
      double PVchi2;
      double PVndof;
      double PVnormalizedChi2;
      double PVx;
      double PVy;
      double PVz;

      // trigger lists
      int _HLT_L1Mu;
      int _HLT_L1MuOpen;
      int _HLT_L2Mu9;
      int _HLT_IsoMu9;
      int _HLT_IsoMu11;
      int _HLT_IsoMu13;
      int _HLT_IsoMu15;
      int _HLT_Mu3;
      int _HLT_Mu5;
      int _HLT_Mu7;
      int _HLT_Mu9;
      int _HLT_Mu11;
      int _HLT_Mu13;
      int _HLT_Mu15;
      int _HLT_Mu15_L1Mu7;
      int _HLT_Mu15_Vtx2cm;   
      int _HLT_Mu15_Vtx2mm;
      int _HLT_DoubleIsoMu3;
      int _HLT_DoubleMu3;
      int _HLT_DoubleMu3_Vtx2cm;
      int _HLT_DoubleMu3_Vtx2mm;
      int _HLT_DoubleMu3_JPsi;
      int _HLT_DoubleMu3_Upsilon;
      int _HLT_DoubleMu7_Z;
      int _HLT_DoubleMu3_SameSign;
      int _HLT_DoubleMu3_Psi2S;


      // Jet
      double JETbDiscriminant[MPSIZE];
      double JETcharge[MPSIZE];
      int JETflavour[MPSIZE];
      int JETntracks[MPSIZE];
      double JETpt[MPSIZE];
      double JETeta[MPSIZE];
      double JETphi[MPSIZE];

      // Pat Muon
      double InvMass[MPSIZE];
      double Angle[MPSIZE];
      double DeltaR[MPSIZE];
      double Dalpha[MPSIZE];
      int isOppSign[MPSIZE];
      int MuonPairIndex[MPSIZE];
      int isHMassPair[MPSIZE];

      double vtxChi2[MPSIZE];
      double vtxPositionX[MPSIZE];
      double vtxPositionY[MPSIZE];
      double vtxPositionZ[MPSIZE];
      double vtxXerror[MPSIZE];
      double vtxYerror[MPSIZE];
      double vtxZerror[MPSIZE];

      int Muon1_muonType[MPSIZE];
      int Muon1_nTrig[MPSIZE];
      int Muon1_triggerObjectType[MPSIZE];
      int Muon1_filterName[MPSIZE];
      double Muon1_phi[MPSIZE];
      double Muon1_eta[MPSIZE];
      double Muon1_pT[MPSIZE];
      double Muon1_Px[MPSIZE];
      double Muon1_Py[MPSIZE];
      double Muon1_Pz[MPSIZE];
      double Muon1_sumtrkpt[MPSIZE];
      double Muon1_trkiso[MPSIZE];
      double Muon1_hcaliso[MPSIZE];
      double Muon1_ecaliso[MPSIZE];
      int Muon1_charge[MPSIZE];
      int Muon1_nChambers[MPSIZE];
      int Muon1_nMatches[MPSIZE];
      int Muon1_stationMask[MPSIZE];
      int Muon1_nSegments[MPSIZE];
      double Muon1_chi2dof[MPSIZE];
      int Muon1_nhits[MPSIZE];
      double Muon1_qoverp[MPSIZE];
      double Muon1_theta[MPSIZE];
      double Muon1_lambda[MPSIZE];
      double Muon1_dxy[MPSIZE];
      double Muon1_d0[MPSIZE];
      double Muon1_dsz[MPSIZE];
      double Muon1_dz[MPSIZE];
      double Muon1_dxyBS[MPSIZE];
      double Muon1_dzBS[MPSIZE];
      double Muon1_dszBS[MPSIZE];
      double Muon1_vx[MPSIZE];
      double Muon1_vy[MPSIZE];
      double Muon1_vz[MPSIZE];
      double Muon1_MCtruth_pT[MPSIZE];
      double Muon1_MCtruth_Px[MPSIZE];
      double Muon1_MCtruth_Py[MPSIZE];
      double Muon1_MCtruth_Pz[MPSIZE];
      double Muon1_MCtruth_eta[MPSIZE];
      double Muon1_MCtruth_phi[MPSIZE];
      int Muon1_MCtruth_charge[MPSIZE];
      int Muon1_MCtruth_mother[MPSIZE];

      int Muon2_muonType[MPSIZE];
      int Muon2_nTrig[MPSIZE];
      int Muon2_triggerObjectType[MPSIZE];
      int Muon2_filterName[MPSIZE];
      double Muon2_phi[MPSIZE];
      double Muon2_eta[MPSIZE];
      double Muon2_pT[MPSIZE];
      double Muon2_Px[MPSIZE];
      double Muon2_Py[MPSIZE];
      double Muon2_Pz[MPSIZE];
      double Muon2_sumtrkpt[MPSIZE];
      double Muon2_trkiso[MPSIZE];
      double Muon2_hcaliso[MPSIZE];
      double Muon2_ecaliso[MPSIZE];
      int Muon2_charge[MPSIZE];
      int Muon2_nChambers[MPSIZE];
      int Muon2_nMatches[MPSIZE];
      int Muon2_stationMask[MPSIZE];
      int Muon2_nSegments[MPSIZE];
      double Muon2_chi2dof[MPSIZE];
      int Muon2_nhits[MPSIZE];
      double Muon2_qoverp[MPSIZE];
      double Muon2_theta[MPSIZE];
      double Muon2_lambda[MPSIZE];
      double Muon2_dxy[MPSIZE];
      double Muon2_d0[MPSIZE];
      double Muon2_dsz[MPSIZE];
      double Muon2_dz[MPSIZE];
      double Muon2_dxyBS[MPSIZE];
      double Muon2_dzBS[MPSIZE];
      double Muon2_dszBS[MPSIZE];
      double Muon2_vx[MPSIZE];
      double Muon2_vy[MPSIZE];
      double Muon2_vz[MPSIZE];
      double Muon2_MCtruth_pT[MPSIZE];
      double Muon2_MCtruth_Px[MPSIZE];
      double Muon2_MCtruth_Py[MPSIZE];
      double Muon2_MCtruth_Pz[MPSIZE];
      double Muon2_MCtruth_eta[MPSIZE];
      double Muon2_MCtruth_phi[MPSIZE];
      int Muon2_MCtruth_charge[MPSIZE];
      int Muon2_MCtruth_mother[MPSIZE];

      // Gen Muon
      double GENInvMass[MPSIZE];
      double GENAngle[MPSIZE];
      int GENisOppSign[MPSIZE];

      double GENMuon1_phi[MPSIZE];
      double GENMuon1_eta[MPSIZE];
      double GENMuon1_pT[MPSIZE];
      double GENMuon1_Px[MPSIZE];
      double GENMuon1_Py[MPSIZE];
      double GENMuon1_Pz[MPSIZE];
      double GENMuon1_mother[MPSIZE];
      int GENMuon1_charge[MPSIZE];

      double GENMuon2_phi[MPSIZE];
      double GENMuon2_eta[MPSIZE];
      double GENMuon2_pT[MPSIZE];
      double GENMuon2_Px[MPSIZE];
      double GENMuon2_Py[MPSIZE];
      double GENMuon2_Pz[MPSIZE];
      double GENMuon2_mother[MPSIZE];
      int GENMuon2_charge[MPSIZE];

};
#endif