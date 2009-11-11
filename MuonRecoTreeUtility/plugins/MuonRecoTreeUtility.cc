
// -*- C++ -*-
//
// Package:    MuonRecoTreeUtility
// Class:      MuonRecoTreeUtility
// 
/**\class MuonRecoTreeUtility MuonRecoTreeUtility.cc Workspace/MuonRecoTreeUtility/src/MuonRecoTreeUtility.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  "Thomas Danielson"
//         Created:  Thu May  8 12:05:03 CDT 2008
// $Id: MuonRecoTreeUtility.cc,v 1.8 2009/11/10 07:25:00 aeverett Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <string>
#include <cstdlib>
#include <utility>
#include <vector>
#include <map>
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackCandidate/interface/TrackCandidate.h>
#include <DataFormats/TrackCandidate/interface/TrackCandidateCollection.h>
#include <DataFormats/TrackingRecHit/interface/TrackingRecHit.h>
#include <DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h>
#include "DataFormats/TrajectoryState/interface/TrackCharge.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"

#include "DataFormats/HLTReco/interface/ModuleTiming.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "TrackingTools/GeomPropagators/interface/StateOnTrackerBound.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

#include "RecoMuon/TrackingTools/interface/MuonErrorMatrix.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include <TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h>

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByPosition.h"

#include <MagneticField/Records/interface/IdealMagneticFieldRecord.h>
#include <MagneticField/Engine/interface/MagneticField.h>

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "TLorentzVector.h"

#include "CommonTools/RecoAlgos/interface/TrackingParticleSelector.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"

// Higher-level muons
#include <DataFormats/TrackReco/interface/Track.h>
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeedCollection.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// so that we can use the pt sorting method, we use this:
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "HepPDT/ParticleDataTable.hh"

#include "RecoMuon/MuonXRay/interface/IDconverttoBinNum.h"
#include "RecoMuon/MuonXRay/interface/MotherSearch.h"
#include "RecoMuon/MuonXRay/interface/DQMHelper.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

//inclusion of muon isolation quantities as seen by the HLT.
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractor.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractorFactory.h"
#include "RecoMuon/L3MuonIsolationProducer/src/L3NominalEfficiencyConfigurator.h"
#include "RecoMuon/MuonIsolation/interface/Cuts.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "RecoMuon/MuonIsolation/interface/Range.h"
//#include "DataFormats/MuonReco/interface/MuonIsolation.h"

#include "RecoMuon/GlobalTrackingTools/interface/GlobalMuonRefitter.h"

#ifdef __CINT__ 
#pragma link C++ class std::map<int,std::vector<double> >++;
#pragma link C++ class std::map<std::string,double>++;
#endif

//
// class decleration
//

class MuonRecoTreeUtility : public edm::EDAnalyzer {
public:
  explicit MuonRecoTreeUtility(const edm::ParameterSet&);
  ~MuonRecoTreeUtility();
  bool isPrimaryMuon(unsigned int inType_)     const { return  inType_ & primaryMuon; }
  bool isSiliconMuon(unsigned int inType_)    const { return  inType_ & siliconMuon; }
  bool isCalConversionMuon(unsigned int inType_) const { return  inType_ & calConversionMuon; }
  bool isOtherMuon(unsigned int inType_) const { return  inType_ & otherMuon; }
  bool isPunchThrough(unsigned int inType_) const { return  inType_ & punchThrough; }
  
 
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
 
  unsigned int getBit(const TrackingParticleRef&) const;
  // ----------member data ---------------------------
  TFile *theFile; // self-explanatory
  TTree *MuTrigData; // reco only.  What you'd fill for actual data.
  TTree *MuTrigMC; // same, but also includes associations, sim muon information.

  //
  // constants, enums and typedefs
  //

  // Members of both MuTrigData and MuTrigMC
  // Event-level information
  int RunNumber;
  int EventNumber;  
  /*
  // Overall MuonHLT execution time, and execution time for each module.
  double totalMuonHLTTime;
  std::map<std::string,double> *muonDigiModuleTimes;
  std::map<std::string,double> *muonLocalRecModuleTimes;
  std::map<std::string,double> *muonL2RecModuleTimes;
  std::map<std::string,double> *muonL3RecModuleTimes;
  std::map<std::string,double> *muonL2IsoModuleTimes;
  std::map<std::string,double> *muonL3IsoModuleTimes;
  std::map<std::string,double> *trackerDigiModuleTimes;
  std::map<std::string,double> *trackerRecModuleTimes;
  std::map<std::string,double> *caloDigiModuleTimes;
  std::map<std::string,double> *caloRecModuleTimes;
  */
  /*
  // Trigger information for the single muon iso and non-iso paths
  // Until they aren't, the level 1 are identical for both...
  int l1SingleMuNonIsoTriggered;
  int l2SingleMuNonIsoTriggered;
  int l3SingleMuNonIsoTriggered;
  int l1SingleMuIsoTriggered;
  int l2SingleMuIsoPreTriggered;
  int l2SingleMuIsoTriggered;
  int l3SingleMuIsoPreTriggered;
  int l3SingleMuIsoTriggered;
  int l1DiMuNonIsoTriggered;
  int l2DiMuNonIsoTriggered;
  int l3DiMuNonIsoTriggered;
  int l1DiMuIsoTriggered;
  int l2DiMuIsoPreTriggered;
  int l2DiMuIsoTriggered;
  int l3DiMuIsoPreTriggered;
  int l3DiMuIsoTriggered;
  std::vector<int> *triggerDecisions;
  std::vector<std::string> *triggerNames;
  */
  // HLT and L1 muon information
  int nMu;
  int nL2;
  int nL3;
  int nTkTracks;
  int nL3Cands;
  //int nL3Seeds;
  // basic L3 muon quantities
  std::vector<int> *muAllGlobalMuons;
  std::vector<int> *muAllStandAloneMuons ;
  std::vector<int> *muAllTrackerMuons ;
  std::vector<int> *muTrackerMuonArbitrated ;
  std::vector<int> *muAllArbitrated ;
  std::vector<int> *muGlobalMuonPromptTight ;
  std::vector<int> *muTMLastStationLoose ;
  std::vector<int> *muTMLastStationTight ;
  std::vector<int> *muTM2DCompatibilityLoose ;
  std::vector<int> *muTM2DCompatibilityTight ;
  std::vector<int> *muTMOneStationLoose ;
  std::vector<int> *muTMOneStationTight ;
  std::vector<int> *muTMLastStationOptimizedLowPtLoose ;
  std::vector<int> *muTMLastStationOptimizedLowPtTight ;
  std::vector<int> *muGMTkChiCompatibility ;
  std::vector<int> *muGMStaChiCompatibility ;
  std::vector<int> *muGMTkKinkTight ;
  std::vector<float> *muCaloCompatibility ;
  std::vector<float> *muSegmentCompatibility ;
  std::vector<float> *muTrkKink ;
  std::vector<float> *muGlbKink ;
  std::vector<float> *muTrkRelChi2 ;
  std::vector<float> *muStaRelChi2 ;
  //  std::vector<int> *muIso03Valid;
  //  std::vector<float> *muIso03sumPt;
  //  std::vector<float> *muIso03emEt;
  //  std::vector<float> *muIso03hadEt;
  //  std::vector<int> *muIso03nTracks;
  //  std::vector<float> *muIso03trackerVetoPt;
  std::vector<int> *muNumberOfChambers;
  std::vector<int> *muNumberOfMatches;
  std::vector<unsigned int> *muStationMask;
  std::map<int,std::vector<int> > *muNCSCSeg;
  std::map<int,std::vector<int> > *muNDTSeg;
  std::map<int,std::vector<int> > *muNRPCSeg;
 

  std::vector<double> *l3P;
  std::vector<double> *l3Px;
  std::vector<double> *l3Py;
  std::vector<double> *l3Pz;
  std::vector<double> *l3Pt;
  std::vector<double> *l3PtError;
  std::vector<double> *l3Pt90;
  std::vector<double> *l3Eta;
  std::vector<double> *l3EtaError;
  std::vector<double> *l3Phi;
  std::vector<double> *l3PhiError;
  std::vector<double> *l3D0;
  std::vector<double> *l3D0Error;
  std::vector<int> *l3NHits;
  std::vector<double> *l3Charge;
  std::vector<double> *l3Chi2;
  std::vector<double> *l3Ndof;
  std::map<int,std::vector<int> > *l3DetIds;
  std::map<int,std::vector<int> > *l3SubdetIds;
  std::map<int,std::vector<int> > *l3Component;  
  std::map<int,std::vector<int> > *l3RecHitsStatus;
  std::map<int,std::vector<double> > *l3RecHitsX;
  std::map<int,std::vector<double> > *l3RecHitsY;
  std::map<int,std::vector<double> > *l3RecHitsZ;
  std::map<int,std::vector<double> > *l3RecHitsXTM;
  std::map<int,std::vector<double> > *l3RecHitsYTM;
  std::map<int,std::vector<double> > *l3RecHitsZTM;
  std::map<int,std::vector<double> > *l3RecHitsXTSOS;
  std::map<int,std::vector<double> > *l3RecHitsYTSOS;
  std::map<int,std::vector<double> > *l3RecHitsZTSOS;
  std::map<int,std::vector<double> > *l3RecHitsPhiTM;
  std::map<int,std::vector<double> > *l3RecHitsErrorTM;
  std::map<int,std::vector<double> > *l3RecHitsPhiTSOS;
  std::map<int, int> *l3NMuHits;
  std::map<int,std::vector<int> > *l3MuStationNumber;
  // L3 muon isolation quantities
  std::vector<double> *l3CalIsoDeposit;
  std::vector<double> *l3TrackIsoDeposit;
  std::vector<double> *l3IsoTrackDR;
  std::vector<double> *l3IsoTrackDRMinDelPt;
  // L3 track fitting and error matrix: note that phi is already in there...
  std::vector<double> *l3Dsz;
  std::vector<double> *l3DszError;
  std::vector<double> *l3Dxy;
  std::vector<double> *l3DxyError;
  std::vector<double> *l3Lambda;
  std::vector<double> *l3LambdaError;
  std::vector<double> *l3Qoverp;
  std::vector<double> *l3QoverpError;
  std::vector<reco::TrackBase::CovarianceMatrix> *l3ErrorMatrix;

  // L2 <-> L3 interface

  // L3 seeding from L2
  std::vector<int> *indexL2SeedingL3;
  std::vector<int> *indexL3SeededFromL2;
  std::vector<int> *l2SeedsL3;

  // Error matrix that's rescaled for the OI algos
  std::map<int,std::vector<double> > *muonErrorMatrix;

  // L3 tracks within window described by muonErrorMatrix

  std::vector<double> *tkTrackP;
  std::vector<double> *tkTrackPx;
  std::vector<double> *tkTrackPy;
  std::vector<double> *tkTrackPz;
  std::vector<double> *tkTrackPt;
  std::vector<double> *tkTrackPtError;
  std::vector<double> *tkTrackEta;
  std::vector<double> *tkTrackEtaError;
  std::vector<double> *tkTrackPhi;
  std::vector<double> *tkTrackPhiError;
  std::vector<double> *tkTrackD0;
  std::vector<double> *tkTrackD0Error;
  std::vector<int> *tkTrackNHits;
  std::vector<double> *tkTrackCharge;
  std::vector<double> *tkTrackChi2;
  std::vector<double> *tkTrackNdof;
  std::map<int,std::vector<int> > *tkTrackDetIds;
  std::map<int,std::vector<int> > *tkTrackSubdetIds;
  std::map<int,std::vector<int> > *tkTrackRecHitsStatus;
  std::map<int,std::vector<double> > *tkTrackRecHitsX;
  std::map<int,std::vector<double> > *tkTrackRecHitsY;
  std::map<int,std::vector<double> > *tkTrackRecHitsZ;
  // L3 track fitting and error matrix: note that phi is already in there...
  std::vector<double> *tkTrackDsz;
  std::vector<double> *tkTrackDszError;
  std::vector<double> *tkTrackDxy;
  std::vector<double> *tkTrackDxyError;
  std::vector<double> *tkTrackLambda;
  std::vector<double> *tkTrackLambdaError;
  std::vector<double> *tkTrackQoverp;
  std::vector<double> *tkTrackQoverpError;
  std::vector<reco::TrackBase::CovarianceMatrix> *tkTrackErrorMatrix;
  
  // basic L2 muon quantities
  std::vector<double> *l2P;
  std::vector<double> *l2Px;
  std::vector<double> *l2Py;
  std::vector<double> *l2Pz;
  std::vector<double> *l2Pt;
  std::vector<double> *l2PtError;
  std::vector<double> *l2Pt90;
  std::vector<double> *l2Eta;
  std::vector<double> *l2EtaError;
  std::vector<double> *l2Phi;
  std::vector<double> *l2PhiError;
  std::vector<double> *l2D0;
  std::vector<double> *l2D0Error;
  std::vector<int> *l2NHits;
  std::vector<double> *l2Charge;
  std::vector<double> *l2Chi2;
  std::vector<double> *l2Ndof;
  std::map<int,std::vector<int> > *l2DetIds;
  std::map<int,std::vector<int> > *l2SubdetIds;
  std::map<int,std::vector<int> > *l2Component;
  std::map<int,std::vector<int> > *l2RecHitsStatus;
  std::map<int,std::vector<double> > *l2RecHitsX;
  std::map<int,std::vector<double> > *l2RecHitsY;
  std::map<int,std::vector<double> > *l2RecHitsZ;
  std::map<int, int> *l2NMuHits;
  std::map<int,std::vector<int> > *l2MuStationNumber;
  // L2 muon isolation quantities
  std::vector<double> *l2CalIsoDeposit;
  // L2 track fitting and error matrix: note that phi is already in there...
  std::vector<double> *l2Dsz;
  std::vector<double> *l2DszError;
  std::vector<double> *l2Dxy;
  std::vector<double> *l2DxyError;
  std::vector<double> *l2Lambda;
  std::vector<double> *l2LambdaError;
  std::vector<double> *l2Qoverp;
  std::vector<double> *l2QoverpError;
  std::vector<reco::TrackBase::CovarianceMatrix> *l2ErrorMatrix;

  // ADAM: associators
  // EXCLUSIVE TO MuTrigMC
  std::vector<int> *l3IsAssociated;
  std::vector<int> *l3ParentID;
  std::vector<int> *l3MotherBinNumber;
  std::vector<double> *l3AssociationVar;
  std::vector<int> *l3AssociationPdgId;
  std::vector<unsigned int> *l3AssociationMyBit;
  std::vector<double> *l3AssociationVtxX;
  std::vector<double> *l3AssociationVtxY;
  std::vector<double> *l3AssociationVtxZ;
  std::vector<int> *l3AssociatedSimMuonIndex;
  std::vector<double> *l3AssociatedSimMuonPt;
  std::vector<double> *l3AssociatedSimMuonEta;
  std::vector<double> *l3AssociatedSimMuonPhi;
  std::vector<int> *l3AssociatedSimMuonNHits;
  std::map<int,std::vector<int> > *l3AssociatedSimMuonDetIds;
  std::map<int,int> *l3AssociatedSimMuonNMuHits;
  std::map<int,std::vector<int> > *l3AssociatedSimMuonMuStationNumber;

  std::vector<double> *l3AssociatedSimMuonDsz;
  std::vector<double> *l3AssociatedSimMuonDxy;
  std::vector<double> *l3AssociatedSimMuonLambda;
  std::vector<double> *l3AssociatedSimMuonQoverP;

  std::vector<int> *tkTrackIsAssociated;
  std::vector<int> *tkTrackParentID;
  std::vector<int> *tkTrackMotherBinNumber;
  std::vector<double> *tkTrackAssociationVar;
  std::vector<int> *tkTrackAssociationPdgId;
  std::vector<unsigned int> *tkTrackAssociationMyBit;
  std::vector<double> *tkTrackAssociationVtxX;
  std::vector<double> *tkTrackAssociationVtxY;
  std::vector<double> *tkTrackAssociationVtxZ;
  std::vector<double> *tkTrackAssociatedSimMuonPt;
  std::vector<double> *tkTrackAssociatedSimMuonEta;
  std::vector<double> *tkTrackAssociatedSimMuonPhi;
  std::vector<int> *tkTrackAssociatedSimMuonNHits;
  std::map<int,std::vector<int> > *tkTrackAssociatedSimMuonDetIds;

  std::vector<double> *tkTrackAssociatedSimMuonDsz;
  std::vector<double> *tkTrackAssociatedSimMuonDxy;
  std::vector<double> *tkTrackAssociatedSimMuonLambda;
  std::vector<double> *tkTrackAssociatedSimMuonQoverP;

  std::vector<int> *l2IsAssociated;
  std::vector<int> *l2ParentID;
  std::vector<int> *l2MotherBinNumber;
  std::vector<double> *l2AssociationVar;
  std::vector<int> *l2AssociationPdgId;
  std::vector<unsigned int> *l2AssociationMyBit;
  std::vector<double> *l2AssociationVtxX;
  std::vector<double> *l2AssociationVtxY;
  std::vector<double> *l2AssociationVtxZ;
  std::vector<double> *l2AssociatedSimMuonPt;
  std::vector<double> *l2AssociatedSimMuonEta;
  std::vector<double> *l2AssociatedSimMuonPhi;
  std::vector<int> *l2AssociatedSimMuonNHits;
  std::map<int,std::vector<int> > *l2AssociatedSimMuonDetIds;
  std::map<int,int> *l2AssociatedSimMuonNMuHits;
  std::map<int,std::vector<int> > *l2AssociatedSimMuonMuStationNumber;

  std::vector<double> *l2AssociatedSimMuonDsz;
  std::vector<double> *l2AssociatedSimMuonDxy;
  std::vector<double> *l2AssociatedSimMuonLambda;
  std::vector<double> *l2AssociatedSimMuonQoverP;

  int nSimMuon;
  std::vector<int> *simMuonParentID;
  std::vector<int> *simMuonMotherBinNumber;
  std::vector<double> *simMuonPt;
  std::vector<double> *simMuonEta;
  std::vector<double> *simMuonPhi;
  std::vector<unsigned int> *simMuonMyBit;
  std::vector<double> *simMuonVtxX;
  std::vector<double> *simMuonVtxY;
  std::vector<double> *simMuonVtxZ;
  std::vector<int> *simMuonNHits;
  std::map<int,std::vector<int> > *simMuonDetIds;
  std::map<int,int> *simMuonNMuHits;
  std::map<int,std::vector<int> > *simMuonMuStationNumber;

  std::vector<double> *simMuonDsz;
  std::vector<double> *simMuonDxy;
  std::vector<double> *simMuonLambda;
  std::vector<double> *simMuonQoverP;

  // SimToReco Associations
  std::vector<int> *simToL3Associated;
  std::vector<double> *simToL3AssociationVar;
  std::vector<int> *simToL3RecoIndex;
  std::vector<int> *simToTkAssociated;
  std::vector<double> *simToTkAssociationVar;
  std::vector<int> *simToTkRecoIndex;
  std::vector<int> *simToL2Associated;
  std::vector<double> *simToL2AssociationVar;
  std::vector<int> *simToL2RecoIndex;

  // And that ends the items going into the tree itself.  Now we just need the things 
  // for putting the items into the tree.
  
  /*
  // Here's where we get the execution time for the modules.
  edm::InputTag theTimerLabel;
  std::vector<std::string> theMuonDigiModules;
  std::vector<std::string> theTrackerDigiModules;
  std::vector<std::string> theTrackerRecModules;
  std::vector<std::string> theTrackerTrackModules;
  std::vector<std::string> theCaloDigiModules;
  std::vector<std::string> theCaloRecModules;
  std::vector<std::string> theMuonLocalRecModules;
  std::vector<std::string> theMuonL2RecModules;
  std::vector<std::string> theMuonL2IsoModules;
  std::vector<std::string> theMuonL3RecModules;
  std::vector<std::string> theMuonL3IsoModules;
  */

  std::string theCategory;
  /*
  std::string singleMuIsoTriggerName;
  std::string singleMuNonIsoTriggerName;
  std::string diMuIsoTriggerName;
  std::string diMuNonIsoTriggerName; 
  */

  std::string outputFileName;

  bool isRecoLevel;

  edm::InputTag l2Label;
  edm::InputTag l3Label;
  edm::InputTag trackLabel;
  edm::InputTag triggerResults_;
  edm::InputTag trackingParticleLabel;
  edm::InputTag allTrackingParticleLabel;
  // To get the seeds for the L2 muons
  edm::InputTag l2SeedCollectionLabel;

  // The all-important error matrix
  edm::ParameterSet errorMatrixPset;
  MuonErrorMatrix * theErrorMatrix;
  std::string thePropagatorName;
  edm::ParameterSet muonServiceParams;
  bool theAdjustAtIp;

  std::string l2AssocLabel;
  std::string l3AssocLabel;
  std::string tkAssocLabel;
  edm::ESHandle<TrackAssociatorBase> l2Associator;
  edm::ESHandle<TrackAssociatorBase> l3Associator;
  edm::ESHandle<TrackAssociatorBase> tkAssociator;

  edm::InputTag theLinkLabel;
  edm::InputTag theMuonLabel;

  edm::ESHandle<MagneticField> field;
  edm::ESHandle<HepPDT::ParticleDataTable> fPDGTable;
  edm::InputTag bsSrc;

  IDconverttoBinNum wantMotherBin;

  //Needed for the ISO quantities
  reco::isodeposit::IsoDepositExtractor* caloDepositExtractor;
  reco::isodeposit::IsoDepositExtractor* trackDepositExtractor;
  reco::isodeposit::IsoDepositExtractor* jetDepositExtractor;
  edm::ParameterSet caloExtractorPSet;
  edm::ParameterSet trackExtractorPSet;
  edm::ParameterSet jetExtractorPSet;
  edm::ParameterSet trackCutsPSet;
  edm::ParameterSet caloCutsPSet;
  edm::InputTag CaloTowerCollectionLabel;
  // This track collection is the tracks for the L3iso track cut
  muonisolation::Cuts L2IsoCalCuts;
  muonisolation::Cuts L3IsoTrackCuts;

  //Needed for the TrackingParticleSelectors
  TrackingParticleSelector tpSelector_primary;
  TrackingParticleSelector tpSelector_silicon;
  TrackingParticleSelector tpSelector_calConversion;

  MuonServiceProxy* theService;
  GlobalMuonRefitter* theGlbRefitter;


public:
  static const unsigned int noBit             =  1<<0;
  static const unsigned int primaryMuon       =  1<<1;
  static const unsigned int siliconMuon       =  1<<2;
  static const unsigned int calConversionMuon =  1<<3;
  static const unsigned int otherMuon         =  1<<4;
  static const unsigned int punchThrough      =  1<<5;


};

//
// static data member definitions
//

double pt90(const reco::TrackRef & tk, const edm::Event & ev){
  std::string prov=ev.getProvenance(tk.id()).moduleLabel();
  std::string instance = ev.getProvenance(tk.id()).productInstanceName();

  double nsigma_Pt=0;

  if (prov.find("standAlone")!= std::string::npos){ nsigma_Pt=3.9;}
  else if (prov.find("global")!= std::string::npos){
    if (instance=="") nsigma_Pt=2.2;
    else nsigma_Pt=1.64; //90% for sure
  }
  else { edm::LogError("p90")<<"provenance: "<<prov<<" instance not recognized.";}

  double pt = tk->pt();
  double err0 = tk->error(0);
  double abspar0 = fabs(tk->parameter(0));
  double ptLx = pt;
  if (abspar0>0) ptLx += nsigma_Pt*err0/abspar0*pt;

  return ptLx;
}

//
// constructors and destructor
//
MuonRecoTreeUtility::MuonRecoTreeUtility(const edm::ParameterSet& iConfig):
  wantMotherBin(iConfig.getParameter<edm::ParameterSet>("IDconverttoBinNum"))
{

  edm::LogInfo("MuonRecoTreeUtility") << "into the constructor.";

  isRecoLevel = iConfig.getUntrackedParameter<bool>("isRecoLevel",false);

  outputFileName = iConfig.getUntrackedParameter<std::string>("outputFileName");

  l2Label = iConfig.getParameter<edm::InputTag>("l2MuonLabel");
  l3Label = iConfig.getParameter<edm::InputTag>("l3MuonLabel");
  trackLabel = iConfig.getParameter<edm::InputTag>("trackLabel");

  triggerResults_ = iConfig.getParameter<edm::InputTag>("triggerResults_");
  trackingParticleLabel = iConfig.getParameter<edm::InputTag>("trackingParticleLabel");
  allTrackingParticleLabel = iConfig.getParameter<edm::InputTag>("allTrackingParticleLabel");
  bsSrc = iConfig.getParameter<edm::InputTag>("beamSpotLabel");

  edm::LogInfo("MuonRecoTreeUtility") << "got the first block of things.  Now getting errorMatrixPset.";
  errorMatrixPset = iConfig.getParameter<edm::ParameterSet>("matrixPset");
  edm::LogInfo("MuonRecoTreeUtility") << "Initializing the matrix rescaler.";
  theErrorMatrix = new MuonErrorMatrix(errorMatrixPset);
  edm::LogInfo("MuonRecoTreeUtility") << "Getting muonServiceParameters.";
  muonServiceParams = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  edm::LogInfo("MuonRecoTreeUtility") << "Getting thePropagatorName.";
  thePropagatorName = iConfig.getParameter<std::string>("propagatorName");
  edm::LogInfo("MuonRecoTreeUtility") << "Getting theAdjustAtIP.";
  theAdjustAtIp = errorMatrixPset.getParameter<bool>("atIP");
  edm::LogInfo("MuonRecoTreeUtility") << "And onwards...";

  theService = new MuonServiceProxy(muonServiceParams);     
      
  // TrackRefitter parameters
  edm::ParameterSet refitterParameters = iConfig.getParameter<edm::ParameterSet>("RefitterParameters");     
  theGlbRefitter = new GlobalMuonRefitter(refitterParameters, theService);
  

  // Get our module timing things
  /*
  theTimerLabel=iConfig.getParameter<edm::InputTag>("TimerLabel"); 
  if(isRecoLevel==false) { LogDebug("MuonRecoTreeUtility");
    edm::ParameterSet ModulesForTiming =  iConfig.getUntrackedParameter<edm::ParameterSet>("TimingModules");
    theMuonDigiModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("MuonDigiModules");
    theMuonLocalRecModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("MuonLocalRecModules");
    theMuonL2RecModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("MuonL2RecModules");
    theMuonL3RecModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("MuonL3RecModules");
    theMuonL2IsoModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("MuonL2IsoModules");
    theMuonL3IsoModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("MuonL3IsoModules");
    theTrackerDigiModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("TrackerDigiModules");
    theTrackerRecModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("TrackerRecModules");
    theCaloDigiModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("CaloDigiModules");
    theCaloRecModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("CaloRecModules");
  }
  */
  /*
  // Get our trigger names from the config
  LogDebug("MuonRecoTreeUtility");
  singleMuIsoTriggerName = iConfig.getParameter<std::string>("singleMuIsoTriggerName");
  singleMuNonIsoTriggerName = iConfig.getParameter<std::string>("singleMuNonIsoTriggerName");
  diMuIsoTriggerName = iConfig.getParameter<std::string>("diMuIsoTriggerName");
  diMuNonIsoTriggerName = iConfig.getParameter<std::string>("diMuNonIsoTriggerName"); 
  */
  LogDebug("MuonRecoTreeUtility");
  l2AssocLabel = iConfig.getParameter<std::string>("l2AssociatorName");
  l3AssocLabel = iConfig.getParameter<std::string>("l3AssociatorName");
  tkAssocLabel = iConfig.getParameter<std::string>("tkAssociatorName");

  theLinkLabel = iConfig.getParameter<edm::InputTag>("linkLabel");
  theMuonLabel = iConfig.getParameter<edm::InputTag>("muonLabel");
  LogDebug("MuonRecoTreeUtility");

  //Make the TrackingParticleSelectors
  edm::ParameterSet tpset_primary = iConfig.getParameter<edm::ParameterSet>("tpSelector_primary");
  tpSelector_primary = 
    TrackingParticleSelector(
			     tpset_primary.getParameter<double>("ptMin"),
			     tpset_primary.getParameter<double>("minRapidity"),
			     tpset_primary.getParameter<double>("maxRapidity"),
			     tpset_primary.getParameter<double>("tip"),
			     tpset_primary.getParameter<double>("lip"),
			     tpset_primary.getParameter<int>("minHit"),
			     tpset_primary.getParameter<bool>("signalOnly"),
			     tpset_primary.getParameter<bool>("chargedOnly"),
			     tpset_primary.getParameter<std::vector<int> >("pdgId"));
  
  edm::ParameterSet tpset_silicon = iConfig.getParameter<edm::ParameterSet>("tpSelector_silicon");
  tpSelector_silicon = 
    TrackingParticleSelector(tpset_silicon.getParameter<double>("ptMin"),
			     tpset_silicon.getParameter<double>("minRapidity"),
			     tpset_silicon.getParameter<double>("maxRapidity"),
			     tpset_silicon.getParameter<double>("tip"),
			     tpset_silicon.getParameter<double>("lip"),
			     tpset_silicon.getParameter<int>("minHit"), 
			     tpset_silicon.getParameter<bool>("signalOnly"),
			     tpset_silicon.getParameter<bool>("chargedOnly"),
			     tpset_silicon.getParameter<std::vector<int> >("pdgId"));
  
  edm::ParameterSet tpset_calConversion = iConfig.getParameter<edm::ParameterSet>("tpSelector_calConversion");
  tpSelector_calConversion = 
    TrackingParticleSelector(tpset_calConversion.getParameter<double>("ptMin"),
			     tpset_calConversion.getParameter<double>("minRapidity"),
			     tpset_calConversion.getParameter<double>("maxRapidity"),
			     tpset_calConversion.getParameter<double>("tip"),
			     tpset_calConversion.getParameter<double>("lip"),
			     tpset_calConversion.getParameter<int>("minHit"),
			     tpset_calConversion.getParameter<bool>("signalOnly"),
			     tpset_calConversion.getParameter<bool>("chargedOnly"),
			     tpset_calConversion.getParameter<std::vector<int> >("pdgId"));
  
  
  // Read in the PSets
  trackCutsPSet = iConfig.getParameter<edm::ParameterSet>("trackCutsPSet");
  caloCutsPSet = iConfig.getParameter<edm::ParameterSet>("caloCutsPSet");

  // Create Extractors to read in the deposits
  LogTrace("MuonRecoTreeUtility")<<"Creating CaloExtractor . . .";
  edm::ParameterSet caloExtractorPSet = iConfig.getParameter<edm::ParameterSet>("CaloExtractorPSet");
  std::string caloExtractorName = caloExtractorPSet.getParameter<std::string>("ComponentName");
  caloDepositExtractor = IsoDepositExtractorFactory::get()->create( caloExtractorName, caloExtractorPSet);

  LogTrace("MuonRecoTreeUtility")<<"Creating TrackExtractor . . .";
  edm::ParameterSet trackExtractorPSet = iConfig.getParameter<edm::ParameterSet>("TrackExtractorPSet");
  std::string trackExtractorName = trackExtractorPSet.getParameter<std::string>("ComponentName");
  trackDepositExtractor = IsoDepositExtractorFactory::get()->create( trackExtractorName, trackExtractorPSet);

  LogTrace("MuonRecoTreeUtility")<<"Creating JetExtractor";
  edm::ParameterSet jetExtractorPSet = iConfig.getParameter<edm::ParameterSet>("JetExtractorPSet");
  std::string jetExtractorName = jetExtractorPSet.getParameter<std::string>("ComponentName");
  jetDepositExtractor = IsoDepositExtractorFactory::get()->create( jetExtractorName, jetExtractorPSet);

  // start in on the L3 cuts
  LogTrace("MuonRecoTreeUtility")<<"Making the Cuts objects . . . ";
  std::string L3IsoTrackCutsName = trackCutsPSet.getParameter<std::string>("L3IsoTrackCutsName");
  if (L3IsoTrackCutsName == "SimpleCuts") {
    LogTrace("MuonRecoTreeUtility")<<". . . with SimpleCuts";
    L2IsoCalCuts = muonisolation::Cuts(caloCutsPSet);
    L3IsoTrackCuts = muonisolation::Cuts(trackCutsPSet);
  }

  // Initialize things so that we have an address for root to write things to
  LogTrace("MuonRecoTreeUtility")<<"Initializing some variables . . . ";
  /*
  triggerDecisions = 0;
  triggerNames = 0;
  */
  /*
  muonDigiModuleTimes = 0;
  muonLocalRecModuleTimes = 0;
  muonL2RecModuleTimes = 0;
  muonL3RecModuleTimes = 0;
  muonL2IsoModuleTimes = 0;
  muonL3IsoModuleTimes = 0;
  trackerDigiModuleTimes = 0;
  trackerRecModuleTimes = 0;
  caloDigiModuleTimes = 0;
  caloRecModuleTimes = 0;
  */

  muAllGlobalMuons = 0;
  muAllStandAloneMuons  = 0;
  muAllTrackerMuons  = 0;
  muTrackerMuonArbitrated  = 0;
  muAllArbitrated  = 0;
  muGlobalMuonPromptTight  = 0;
  muTMLastStationLoose  = 0;
  muTMLastStationTight  = 0;
  muTM2DCompatibilityLoose  = 0;
  muTM2DCompatibilityTight  = 0;
  muTMOneStationLoose  = 0;
  muTMOneStationTight  = 0;
  muTMLastStationOptimizedLowPtLoose  = 0;
  muTMLastStationOptimizedLowPtTight  = 0;
  muGMTkChiCompatibility  = 0;
  muGMStaChiCompatibility  = 0;
  muGMTkKinkTight  = 0;
  muCaloCompatibility  = 0;
  muSegmentCompatibility  = 0;
  muTrkKink  = 0;
  muGlbKink  = 0;
  muTrkRelChi2  = 0;
  muStaRelChi2 = 0;
  //  muIso03Valid = 0;
  //  muIso03sumPt = 0;
  //  muIso03emEt = 0;
  //  muIso03hadEt = 0;
  //  muIso03nTracks = 0;
  //  muIso03trackerVetoPt = 0;
  muNumberOfChambers =0;
  muNumberOfMatches =0;
  muStationMask =0;
  muNCSCSeg =  new std::map<int,std::vector<int> >;
  muNDTSeg =  new std::map<int,std::vector<int> >;
  muNRPCSeg =  new std::map<int,std::vector<int> >;
  l3P = 0;
  l3Pt = 0;
  l3Px = 0;
  l3Py = 0;
  l3Pz = 0;
  l3PtError = 0;
  l3Pt90 = 0;
  l3Eta = 0;
  l3EtaError = 0;
  l3Phi = 0;
  l3PhiError = 0;
  l3D0 = 0;
  l3D0Error = 0;
  l3NHits = 0;
  l3Charge = 0;
  l3Chi2 = 0;
  l3Ndof = 0;
  l3DetIds = new std::map<int,std::vector<int> >;
  l3SubdetIds = new std::map<int,std::vector<int> >;
  l3NMuHits = new std::map<int,int>;
  l3MuStationNumber = new std::map<int, std::vector <int> >;
  l3Component = new std::map<int,std::vector<int> >;
  l3RecHitsStatus = new std::map<int,std::vector<int> >;
  l3RecHitsX = new std::map<int,std::vector<double> >;
  l3RecHitsY = new std::map<int,std::vector<double> >;
  l3RecHitsZ = new std::map<int,std::vector<double> >;
  l3RecHitsXTM = new std::map<int,std::vector<double> >;
  l3RecHitsYTM = new std::map<int,std::vector<double> >;
  l3RecHitsZTM = new std::map<int,std::vector<double> >;
  l3RecHitsXTSOS = new std::map<int,std::vector<double> >;
  l3RecHitsYTSOS = new std::map<int,std::vector<double> >;
  l3RecHitsZTSOS = new std::map<int,std::vector<double> >;
  l3RecHitsPhiTM = new std::map<int,std::vector<double> >;
  l3RecHitsErrorTM = new std::map<int,std::vector<double> >;
  l3RecHitsPhiTSOS = new std::map<int,std::vector<double> >;
  
  //l3CalIsoDeposit = 0;
  //l3TrackIsoDeposit = 0;
  l3IsoTrackDR = 0;
  l3IsoTrackDRMinDelPt = 0;
  
  l3Dsz = 0;
  l3DszError = 0;
  l3Dxy = 0;
  l3DxyError = 0;
  l3Lambda = 0;
  l3LambdaError = 0;
  l3Qoverp = 0;
  l3QoverpError = 0;
  l3ErrorMatrix = 0;
  
  indexL2SeedingL3 = 0;
  indexL3SeededFromL2 = 0;
  l2SeedsL3 = 0;
  LogDebug("MuonRecoTreeUtility");
  muonErrorMatrix = new std::map<int,std::vector<double> >;

  l3IsAssociated = 0;
  l3ParentID = 0;
  l3MotherBinNumber = 0;
  l3AssociationVar = 0;
  l3AssociationPdgId = 0;
  l3AssociationMyBit = 0;
  l3AssociationVtxX = 0;
  l3AssociationVtxY = 0;
  l3AssociationVtxZ = 0;
  l3AssociatedSimMuonIndex = 0;
  l3AssociatedSimMuonPt = 0;
  l3AssociatedSimMuonEta = 0;
  l3AssociatedSimMuonPhi = 0;
  l3AssociatedSimMuonNHits = 0;
  l3AssociatedSimMuonDetIds = new std::map<int,std::vector<int> >;
  l3AssociatedSimMuonNMuHits = new std::map<int, int>;
  l3AssociatedSimMuonMuStationNumber = new std::map<int, std::vector<int> >;

  l3AssociatedSimMuonDsz = 0;
  l3AssociatedSimMuonDxy = 0;
  l3AssociatedSimMuonLambda = 0;
  l3AssociatedSimMuonQoverP = 0;

  tkTrackP = 0;
  tkTrackPt = 0;
  tkTrackPx = 0;
  tkTrackPy = 0;
  tkTrackPz = 0;
  tkTrackPtError = 0;
  tkTrackEta = 0;
  tkTrackEtaError = 0;
  tkTrackPhi = 0;
  tkTrackPhiError = 0;
  tkTrackD0 = 0;
  tkTrackD0Error = 0;
  tkTrackNHits = 0;
  tkTrackCharge = 0;
  tkTrackChi2 = 0;
  tkTrackNdof = 0;
  tkTrackDetIds = new std::map<int,std::vector<int> >;
  tkTrackSubdetIds = new std::map<int,std::vector<int> >;
  tkTrackRecHitsStatus = new std::map<int,std::vector<int> >;
  tkTrackRecHitsX = new std::map<int,std::vector<double> >;
  tkTrackRecHitsY = new std::map<int,std::vector<double> >;
  tkTrackRecHitsZ = new std::map<int,std::vector<double> >;

  tkTrackDsz = 0;
  tkTrackDszError = 0;
  tkTrackDxy = 0;
  tkTrackDxyError = 0;
  tkTrackLambda = 0;
  tkTrackLambdaError = 0;
  tkTrackQoverp = 0;
  tkTrackQoverpError = 0;
  tkTrackErrorMatrix = 0;

  tkTrackIsAssociated = 0;
  tkTrackParentID = 0;
  tkTrackMotherBinNumber = 0;
  tkTrackAssociationVar = 0;
  tkTrackAssociationPdgId = 0;
  tkTrackAssociationMyBit = 0;
  tkTrackAssociationVtxX = 0;
  tkTrackAssociationVtxY = 0;
  tkTrackAssociationVtxZ = 0;
  tkTrackAssociatedSimMuonPt = 0;
  tkTrackAssociatedSimMuonEta = 0;
  tkTrackAssociatedSimMuonPhi = 0;
  tkTrackAssociatedSimMuonNHits = 0;
  tkTrackAssociatedSimMuonDetIds = new std::map<int,std::vector<int> >;

  tkTrackAssociatedSimMuonDsz = 0;
  tkTrackAssociatedSimMuonDxy = 0;
  tkTrackAssociatedSimMuonLambda = 0;
  tkTrackAssociatedSimMuonQoverP = 0;

  l2P = 0;
  l2Px = 0;
  l2Py = 0;
  l2Pz = 0;
  l2Pt = 0;
  l2PtError = 0;
  l2Pt90 = 0;
  l2Eta = 0;
  l2EtaError = 0;
  l2Phi = 0;
  l2PhiError = 0;
  l2D0 = 0;
  l2D0Error = 0;
  l2NHits = 0;
  l2Charge = 0;
  l2Chi2 = 0;
  l2Ndof = 0;
  l2DetIds = new std::map<int,std::vector<int> >;
  l2SubdetIds = new std::map<int,std::vector<int> >;
  l2Component = new std::map<int,std::vector<int> >;
  l2NMuHits = new std::map<int,int>;
  l2MuStationNumber = new std::map<int, std::vector <int> >;
  l2RecHitsStatus = new std::map<int,std::vector<int> >;
  l2RecHitsX = new std::map<int,std::vector<double> >;
  l2RecHitsY = new std::map<int,std::vector<double> >;
  l2RecHitsZ = new std::map<int,std::vector<double> >;

  l2CalIsoDeposit = 0;

  l2Dsz = 0;
  l2DszError = 0;
  l2Dxy = 0;
  l2DxyError = 0;
  l2Lambda = 0;
  l2LambdaError = 0;
  l2Qoverp = 0;
  l2QoverpError = 0;
  l2ErrorMatrix = 0;

  l2IsAssociated = 0;
  l2ParentID = 0;
  l2MotherBinNumber = 0;
  l2AssociationVar = 0;
  l2AssociationPdgId = 0;
  l2AssociationMyBit = 0;
  l2AssociationVtxX = 0;
  l2AssociationVtxY = 0;
  l2AssociationVtxZ = 0;
  l2AssociatedSimMuonPt = 0;
  l2AssociatedSimMuonEta = 0;
  l2AssociatedSimMuonPhi = 0;
  l2AssociatedSimMuonNHits = 0;
  l2AssociatedSimMuonDetIds = new std::map<int,std::vector<int> >;
  l2AssociatedSimMuonNMuHits = new std::map<int, int>;
  l2AssociatedSimMuonMuStationNumber = new std::map<int, std::vector<int> >;

  l2AssociatedSimMuonDsz = 0;
  l2AssociatedSimMuonDxy = 0;
  l2AssociatedSimMuonLambda = 0;
  l2AssociatedSimMuonQoverP = 0;

  simMuonParentID = 0;
  simMuonMotherBinNumber = 0;
  simMuonPt = 0;
  simMuonEta = 0;
  simMuonPhi = 0;
  simMuonMyBit = 0;
  simMuonVtxX = 0;
  simMuonVtxY = 0;
  simMuonVtxZ = 0;
  simMuonNHits = 0;
  simMuonDetIds = new std::map<int,std::vector<int> >;
  simMuonNMuHits = new std::map<int, int>;
  simMuonMuStationNumber = new std::map<int, std::vector<int> >;
    

  simMuonDsz = 0;
  simMuonDxy = 0;
  simMuonLambda = 0;
  simMuonQoverP = 0;

  simToL3Associated = 0;
  simToL3AssociationVar = 0;
  simToL3RecoIndex = 0;
  simToTkAssociated = 0;
  simToTkAssociationVar = 0;
  simToTkRecoIndex = 0;
  simToL2Associated = 0;
  simToL2AssociationVar = 0;
  simToL2RecoIndex = 0;


}


MuonRecoTreeUtility::~MuonRecoTreeUtility()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  if (theService) delete theService;
  if (theGlbRefitter) delete theGlbRefitter;

}


//
// member functions
//

// ------------ method called to for each event  ------------
void MuonRecoTreeUtility::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;

  edm::LogInfo("MuonRecoTreeUtility") << "Beginning of the loop.";

  theService->update(iSetup);

  theGlbRefitter->setEvent(iEvent);

  theGlbRefitter->setServices(theService->eventSetup());


  //get the mag field and the beamspot
  iSetup.get<IdealMagneticFieldRecord>().get(field);
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel(bsSrc,recoBeamSpotHandle);
  reco::BeamSpot bs = *recoBeamSpotHandle;      

  /*
  edm::Handle<EventTime> evtTime;
  iEvent.getByLabel(theTimerLabel, evtTime); 
  */

  //get the collection of muons at all levels
  edm::Handle<reco::TrackCollection> l2Muons;
  edm::Handle<reco::TrackCollection> l3Muons;
  edm::Handle<reco::TrackCollection> tkTracks;

  iEvent.getByLabel(l2Label,l2Muons);
  iEvent.getByLabel(l3Label,l3Muons);
  iEvent.getByLabel(trackLabel,tkTracks);

  edm::Handle<edm::View<reco::Track> > l2MuonsForAssociation;
  edm::Handle<edm::View<reco::Track> > l3MuonsForAssociation;
  edm::Handle<edm::View<reco::Track> > tkTracksForAssociation;
  iEvent.getByLabel(l2Label,l2MuonsForAssociation);
  iEvent.getByLabel(l3Label,l3MuonsForAssociation);
  iEvent.getByLabel(trackLabel,tkTracksForAssociation);

  //open a collection of Tracking Particles
  edm::Handle<TrackingParticleCollection> TPtracks;
  iEvent.getByLabel(trackingParticleLabel,TPtracks);
  edm::Handle<TrackingParticleCollection> allTPtracks;
  iEvent.getByLabel(allTrackingParticleLabel,allTPtracks);

  //get a hold off simulated information
  Handle<SimTrackContainer> SimTk;
  Handle<SimVertexContainer> SimVtx;
  //get a hold on generated information
  Handle<HepMCProduct> hepmc;

  iEvent.getByType(hepmc);
  iEvent.getByLabel("g4SimHits",SimVtx);
  iEvent.getByLabel("g4SimHits",SimTk);

  //get the associators
  iSetup.get<TrackAssociatorRecord>().get(l2AssocLabel,l2Associator);
  iSetup.get<TrackAssociatorRecord>().get(l3AssocLabel,l3Associator);
  //  iSetup.get<TrackAssociatorRecord>().get(k3AssocLabel,l3Associator);

  //ADAM: associator
  //associate RecoToSim
  reco::RecoToSimCollection l2RecSimColl = l2Associator->associateRecoToSim(l2MuonsForAssociation, TPtracks, &iEvent);
  reco::RecoToSimCollection l3RecSimColl = l3Associator->associateRecoToSim(l3MuonsForAssociation, TPtracks, &iEvent);
  reco::RecoToSimCollection tkRecSimColl = l3Associator->associateRecoToSim(l3MuonsForAssociation, TPtracks, &iEvent);
  reco::RecoToSimCollection l3RecAllSimColl = l3Associator->associateRecoToSim(l3MuonsForAssociation, allTPtracks, &iEvent);

  //associate SimToReco
  reco::SimToRecoCollection l2SimRecColl = l2Associator->associateSimToReco(l2MuonsForAssociation, TPtracks, &iEvent);
  reco::SimToRecoCollection l3SimRecColl = l3Associator->associateSimToReco(l3MuonsForAssociation, TPtracks, &iEvent);
  reco::SimToRecoCollection tkSimRecColl = l3Associator->associateSimToReco(l3MuonsForAssociation, TPtracks, &iEvent);
  reco::SimToRecoCollection l3AllSimRecColl = l3Associator->associateSimToReco(l3MuonsForAssociation, allTPtracks, &iEvent);

  // Event-level information: run, event, and Trigger Table
  EventNumber = iEvent.id().event();
  RunNumber = iEvent.id().run();   
  LogDebug("MuonRecoTreeUtility")<<"Run " << RunNumber << " Event  " << EventNumber;
  // Begin entering muon quantities.  Number of muons at each level.
  // These will be the loop limits over L1, L2, and L3.
  nL2 = l2Muons->size();
  nL3 = l3Muons->size();
  if (!tkTracks.failedToGet()) {
    nTkTracks = tkTracks->size();
  }
  else {
    nTkTracks = -1;
    edm::LogInfo("MuonRecoTreeUtility") << "no l3 tracks";
  }

  //  edm::LogInfo("MuonRecoTreeUtility") << "How many L1, L2, L3 do we have? " << nL1 << " " << nL2 << " " << nL3;

  //ISO variables go here
  LogDebug("MuonRecoTreeUtility");
  caloDepositExtractor->fillVetos(iEvent,iSetup,*l2Muons);
  trackDepositExtractor->fillVetos(iEvent,iSetup,*l3Muons);
  edm::LogInfo("MuonRecoTreeUtility") << "veto filled";

  LogDebug("MuonRecoTreeUtility");
  reco::IsoDeposit::Vetos trackVetos;
  typedef std::vector< std::pair<reco::TrackRef,reco::IsoDeposit> > MuonsWithDeposits;
  MuonsWithDeposits muonsWithDeposits;
  LogDebug("MuonRecoTreeUtility");
  edm::LogInfo("MuonRecoTreeUtility") << "deposit objects created.";

  // get hold of the TriggerResults objects
  /*
  edm::Handle<TriggerResults> triggerResults;
  iEvent.getByLabel(triggerResults_,triggerResults);
  TriggerNames namesOfTriggers(*triggerResults);
  */

  // Right, that's the setup done.  Now let's put in our execution times.
  /*
  edm::LogInfo("MuonRecoTreeUtility") << "doing the module timing.";

  unsigned nTotalModules = evtTime->size();
  totalMuonHLTTime = 0;
  for(unsigned int i = 0; i != nTotalModules; ++i){
    std::string module_name = evtTime->name(i);
    for ( unsigned int j = 0; j != theMuonDigiModules.size(); ++j ) {
      if ( theMuonDigiModules[j] == module_name) {
	totalMuonHLTTime+=evtTime->time(i);
	(*muonDigiModuleTimes).insert(std::make_pair(theMuonDigiModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theMuonLocalRecModules.size(); ++j ) {
      if ( theMuonLocalRecModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
        (*muonLocalRecModuleTimes).insert(std::make_pair(theMuonLocalRecModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theMuonL2RecModules.size(); ++j ) {
      if ( theMuonL2RecModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
	(*muonL2RecModuleTimes).insert(std::make_pair(theMuonL2RecModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theMuonL3RecModules.size(); ++j ) {
      if ( theMuonL3RecModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
        (*muonL3RecModuleTimes).insert(std::make_pair(theMuonL3RecModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theMuonL2IsoModules.size(); ++j ) {
      if ( theMuonL2IsoModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
        (*muonL2IsoModuleTimes).insert(std::make_pair(theMuonL2IsoModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theMuonL3IsoModules.size(); ++j ) {
      if ( theMuonL3IsoModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
        (*muonL3IsoModuleTimes).insert(std::make_pair(theMuonL3IsoModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theTrackerDigiModules.size(); ++j ) {
      if ( theTrackerDigiModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
        (*trackerDigiModuleTimes).insert(std::make_pair(theTrackerDigiModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theTrackerRecModules.size(); ++j ) {
      if ( theTrackerRecModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
        (*trackerRecModuleTimes).insert(std::make_pair(theTrackerRecModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theCaloDigiModules.size(); ++j ) {
      if ( theCaloDigiModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
	(*caloDigiModuleTimes).insert(std::make_pair(theCaloDigiModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theCaloRecModules.size(); ++j ) {
      if ( theCaloRecModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
        (*caloRecModuleTimes).insert(std::make_pair(theCaloRecModules[j],evtTime->time(i)));
      }
    }
  }
  */

  /*
  edm::LogInfo("MuonRecoTreeUtility") << "checking what filter passed.";

  unsigned int indexSingleMuIso = namesOfTriggers.triggerIndex(singleMuIsoTriggerName );
  unsigned int indexSingleMuNoIso = namesOfTriggers.triggerIndex(singleMuNonIsoTriggerName);
  unsigned int indexDiMuIso = namesOfTriggers.triggerIndex(diMuIsoTriggerName );
  unsigned int indexDiMuNoIso = namesOfTriggers.triggerIndex(diMuNonIsoTriggerName); 
  //  edm::LogInfo("MuonRecoTreeUtility") << "Trying to get the trigger decisions.";
  edm::LogInfo("MuonRecoTreeUtility") << "Indexes: " << indexSingleMuIso << " " << indexSingleMuNoIso << " " << indexDiMuIso << " " << indexDiMuNoIso;

  if (indexSingleMuNoIso<namesOfTriggers.size()){
    if (triggerResults->index(indexSingleMuNoIso) > 6) l1SingleMuNonIsoTriggered = 1;
    else l1SingleMuNonIsoTriggered = 0;
    if (triggerResults->index(indexSingleMuNoIso) > 20) l2SingleMuNonIsoTriggered = 1;
    else l2SingleMuNonIsoTriggered = 0;
    if (triggerResults->state(indexSingleMuNoIso) == 1) l3SingleMuNonIsoTriggered = 1;
    else l3SingleMuNonIsoTriggered = 0;
  }
  else
    edm::LogError("MuonRecoTreeUtility")<<singleMuNonIsoTriggerName<<" is not a valid trigger name here.";
  
  if (indexSingleMuIso<namesOfTriggers.size()){
    if (triggerResults->index(indexSingleMuIso) > 6) l1SingleMuIsoTriggered = 1;
    else l1SingleMuIsoTriggered = 0;
    if (triggerResults->index(indexSingleMuIso) > 20) l2SingleMuIsoPreTriggered = 1;
    else l2SingleMuIsoPreTriggered = 0;
    if (triggerResults->index(indexSingleMuIso) > 34) l2SingleMuIsoTriggered = 1;
    else l2SingleMuIsoTriggered = 0;
    if (triggerResults->index(indexSingleMuIso) > 44) l3SingleMuIsoPreTriggered = 1;
    else l3SingleMuIsoPreTriggered = 0;
    if (triggerResults->state(indexSingleMuIso) == 1) l3SingleMuIsoTriggered = 1;
    else l3SingleMuIsoTriggered = 0;
  }
  else
    edm::LogError("MuonRecoTreeUtility")<<singleMuIsoTriggerName<<" is not a valid trigger name here.";
  
  if (indexDiMuIso<namesOfTriggers.size()){
    if (triggerResults->index(indexDiMuNoIso) > 6) l1DiMuNonIsoTriggered = 1;
    else l1DiMuNonIsoTriggered = 0;
    if (triggerResults->index(indexDiMuNoIso) > 20) l2DiMuNonIsoTriggered = 1;
    else l2DiMuNonIsoTriggered = 0;
    if (triggerResults->state(indexDiMuNoIso) == 1) l3DiMuNonIsoTriggered = 1;
    else l3DiMuNonIsoTriggered = 0;
  }
  else
    edm::LogError("MuonRecoTreeUtility")<<diMuNonIsoTriggerName<<" is not a valid trigger name here.";
  
  if (indexDiMuIso<namesOfTriggers.size()){
    if (triggerResults->index(indexDiMuIso) > 6) l1DiMuIsoTriggered = 1;
    else l1DiMuIsoTriggered = 0;
    if (triggerResults->index(indexDiMuIso) > 20) l2DiMuIsoPreTriggered = 1;
    else l2DiMuIsoPreTriggered = 0;
    if (triggerResults->index(indexDiMuIso) > 34) l2DiMuIsoTriggered = 1;
    else l2DiMuIsoTriggered = 0;
    if (triggerResults->index(indexDiMuIso) > 44) l3DiMuIsoPreTriggered = 1;
    else l3DiMuIsoPreTriggered = 0;
    if (triggerResults->state(indexDiMuIso) == 1) l3DiMuIsoTriggered = 1;
    else l3DiMuIsoTriggered = 0;
  }
  else
    edm::LogError("MuonRecoTreeUtility")<<diMuIsoTriggerName<<" is not a valid trigger name here.";
  
  edm::LogInfo("MuonRecoTreeUtility") << "done.";

  for (unsigned int i = 0; i < triggerResults->size(); i++) {
    (*triggerDecisions).push_back(triggerResults->state(i));
    (*triggerNames).push_back(namesOfTriggers.triggerName(i));
  }
  */ 

  //edm::Handle<reco::MuonTrackLinksCollection> l3ToL2Links;
  //iEvent.getByLabel(theLinkLabel, l3ToL2Links);

  //start ADAM recoMuon
  //edm::Handle<reco::MuonCollection> recMuons;
  //iEvent.getByLabel(theMuonLabel, recMuons);

  // Get Muon Tracks
  //Handle<View<Track> > trkMuHandle;
  //iEvent.getByLabel(trkMuLabel_, trkMuHandle);
  //View<Track> trkMuColl = *(trkMuHandle.product());

  //Handle<View<Track> > staMuHandle;
  //iEvent.getByLabel(staMuLabel_, staMuHandle);
  //View<Track> staMuColl = *(staMuHandle.product());

  //Handle<View<Track> > glbMuHandle;
  //iEvent.getByLabel(glbMuLabel_, glbMuHandle);
  //View<Track> glbMuColl = *(glbMuHandle.product());

  // Get Muons
  Handle<View<Muon> > muonHandle;
  iEvent.getByLabel(theMuonLabel, muonHandle);
  View<Muon> muonColl = *(muonHandle.product());

  nMu = muonColl.size();

  int iMu = 0;
  for(View<Muon>::const_iterator iMuon = muonColl.begin();
      iMuon != muonColl.end(); ++iMuon) {

    //start the recoMuon member block
    {
      (*muAllGlobalMuons).push_back( muon::isGoodMuon(*iMuon,muon::AllGlobalMuons) );
      (*muAllStandAloneMuons ).push_back( muon::isGoodMuon(*iMuon,muon::AllStandAloneMuons));
      (*muAllTrackerMuons ).push_back( muon::isGoodMuon(*iMuon,muon::AllTrackerMuons));
      (*muTrackerMuonArbitrated ).push_back( muon::isGoodMuon(*iMuon,muon::TrackerMuonArbitrated));
      (*muAllArbitrated ).push_back( muon::isGoodMuon(*iMuon,muon::AllArbitrated));
      (*muGlobalMuonPromptTight ).push_back( muon::isGoodMuon(*iMuon,muon::GlobalMuonPromptTight));
      (*muTMLastStationLoose ).push_back( muon::isGoodMuon(*iMuon,muon::TMLastStationLoose));
      (*muTMLastStationTight ).push_back( muon::isGoodMuon(*iMuon,muon::TMLastStationTight));
      (*muTM2DCompatibilityLoose ).push_back( muon::isGoodMuon(*iMuon,muon::TM2DCompatibilityLoose));
      (*muTM2DCompatibilityTight ).push_back( muon::isGoodMuon(*iMuon,muon::TM2DCompatibilityTight));
      (*muTMOneStationLoose ).push_back( muon::isGoodMuon(*iMuon,muon::TMOneStationLoose));
      (*muTMOneStationTight ).push_back( muon::isGoodMuon(*iMuon,muon::TMOneStationTight));
      (*muTMLastStationOptimizedLowPtLoose ).push_back( muon::isGoodMuon(*iMuon,muon::TMLastStationOptimizedLowPtLoose));
      (*muTMLastStationOptimizedLowPtTight ).push_back( muon::isGoodMuon(*iMuon,muon::TMLastStationOptimizedLowPtTight));
      (*muGMTkChiCompatibility ).push_back( muon::isGoodMuon(*iMuon,muon::GMTkChiCompatibility));
      (*muGMStaChiCompatibility ).push_back( muon::isGoodMuon(*iMuon,muon::GMStaChiCompatibility));
      (*muGMTkKinkTight ).push_back( muon::isGoodMuon(*iMuon,muon::GMTkKinkTight));

      (*muCaloCompatibility).push_back(muon::caloCompatibility(*iMuon));
      (*muSegmentCompatibility).push_back(muon::segmentCompatibility(*iMuon));

      (*muNumberOfChambers).push_back(iMuon->numberOfChambers());
      (*muNumberOfMatches).push_back(iMuon->numberOfMatches());
      (*muStationMask).push_back(iMuon->stationMask(Muon::SegmentArbitration));

      std::vector<int> *nCSCSeg = new std::vector<int>;
      std::vector<int> *nDTSeg = new std::vector<int>;
      std::vector<int> *nRPCSeg = new std::vector<int>;
      for(int station = 0; station < 4; ++station) {
	nCSCSeg->push_back(iMuon->numberOfSegments(station+1, MuonSubdetId::CSC, Muon::NoArbitration));
	nDTSeg->push_back(iMuon->numberOfSegments(station+1, MuonSubdetId::DT, Muon::NoArbitration));
	nRPCSeg->push_back(iMuon->numberOfSegments(station+1, MuonSubdetId::RPC, Muon::NoArbitration));
      }

      (*muNCSCSeg).insert(std::make_pair(iMu,*nCSCSeg));
      (*muNDTSeg).insert(std::make_pair(iMu,*nDTSeg));
      (*muNRPCSeg).insert(std::make_pair(iMu,*nRPCSeg));

      //      bool isoValid = iMuon->isIsolationValid();
      //      MuonIsolation muIso03 = iMuon->isolationR03();
      //      float sumPt = isoValid ? muIso03.sumPt : -999.;
      //      float emEt = isoValid ? muIso03.emEt : -999.;
      //      float hadEt = isoValid ? muIso03.hadEt : -999.;
      //      int nTracks = isoValid ? muIso03.nTracks : 0;
      //      float trackerVetoPt = isoValid ? muIso03.trackerVetoPt : -999.;
      //      (*muIso03Valid).push_back(isoValid);
      //      (*muIso03sumPt).push_back(sumPt);
      //      (*muIso03emEt).push_back(emEt);
      //      (*muIso03hadEt).push_back(hadEt);
      //      (*muIso03nTracks).push_back(nTracks);
      //      (*muIso03trackerVetoPt).push_back(trackerVetoPt);

      if(iMuon->isGlobalMuon() && iMuon->isQualityValid()) {
	(*muTrkKink).push_back(iMuon->combinedQuality().trkKink);
	(*muGlbKink).push_back(iMuon->combinedQuality().glbKink);
	(*muTrkRelChi2).push_back(iMuon->combinedQuality().trkRelChi2);
	(*muStaRelChi2).push_back(iMuon->combinedQuality().staRelChi2);
      } else {
	(*muTrkKink).push_back(-999.);
	(*muGlbKink).push_back(-999.);
	(*muTrkRelChi2).push_back(-999.);
	(*muStaRelChi2).push_back(-999.);
      }
    }

    const reco::TrackRef glbTrack = ( iMuon->isGlobalMuon()) ? 
      iMuon->combinedMuon() : reco::TrackRef();

    const RefToBase<Track> glbTrackRB(glbTrack);
    
    const reco::TrackRef tkTrack = ( iMuon->isTrackerMuon() ) ? 
      iMuon->innerTrack() : TrackRef();
    
    const reco::TrackRef staTrack = ( iMuon->isStandAloneMuon() ) ? 
      iMuon->outerTrack() : TrackRef();

    if(glbTrack.isAvailable()) {
      int iL3 = iMu;
      const reco::TrackRef refL3(glbTrack);

      // Fill in basic information of l3Muons
      (*l3P).push_back(refL3->p());
      (*l3Px).push_back(refL3->px());
      (*l3Py).push_back(refL3->py());
      (*l3Pz).push_back(refL3->pz());
      (*l3Pt).push_back(refL3->pt());
      (*l3PtError).push_back(refL3->ptError());
      (*l3Pt90).push_back(pt90(refL3,iEvent));
      (*l3Eta).push_back(refL3->eta());
      (*l3EtaError).push_back(refL3->etaError());
      (*l3Phi).push_back(refL3->phi());
      (*l3PhiError).push_back(refL3->phiError());
      (*l3D0).push_back(refL3->d0());
      (*l3D0Error).push_back(refL3->d0Error());
      (*l3NHits).push_back(refL3->recHitsSize());
      (*l3Charge).push_back(refL3->charge());
      (*l3Chi2).push_back(refL3->chi2());
      (*l3Ndof).push_back(refL3->ndof());
      // Fill in the track fitting parameters (with phi filled in above)
      (*l3Dsz).push_back(refL3->dsz());
      (*l3DszError).push_back(refL3->dszError());
      (*l3Dxy).push_back(refL3->dxy());
      (*l3DxyError).push_back(refL3->dxyError());
      (*l3Lambda).push_back(refL3->lambda());
      (*l3LambdaError).push_back(refL3->lambdaError());
      (*l3Qoverp).push_back(refL3->qoverp());
      (*l3QoverpError).push_back(refL3->qoverpError());
      (*l3ErrorMatrix).push_back(refL3->covariance());

      std::vector<int> *idsForThisL3 = new std::vector<int>;
      std::vector<int> *subidsForThisL3 = new std::vector<int>;
      std::vector<int> *detsForThisL3 = new std::vector<int>;
      std::vector<int> *statusForThisL3 = new std::vector<int>;
      std::vector<double> *xForThisL3 = new std::vector<double>;
      std::vector<double> *yForThisL3 = new std::vector<double>;
      std::vector<double> *zForThisL3 = new std::vector<double>;
      std::vector<double> *xForThisL3TM = new std::vector<double>;
      std::vector<double> *yForThisL3TM = new std::vector<double>;
      std::vector<double> *zForThisL3TM = new std::vector<double>;
      std::vector<double> *xForThisL3TSOS = new std::vector<double>;
      std::vector<double> *yForThisL3TSOS = new std::vector<double>;
      std::vector<double> *zForThisL3TSOS = new std::vector<double>;
      std::vector<double> *errorForThisL3TM = new std::vector<double>;
      std::vector<double> *phiForThisL3TM = new std::vector<double>;
      std::vector<double> *phiForThisL3TSOS = new std::vector<double>;
      std::vector<int> *stationsForThisL3 = new std::vector<int>;
      int nMuHitsForThisL3 = 0;

      edm::ESHandle<TransientTrackingRecHitBuilder> trackBuilder;
      edm::ESHandle<TransientTrackingRecHitBuilder> muonBuilder;
      std::string trackBuilderName = "WithTrackAngle";
      std::string muonBuilderName = "MuonRecHitBuilder";
      iSetup.get<TransientRecHitRecord>().get(trackBuilderName,trackBuilder);
      iSetup.get<TransientRecHitRecord>().get(muonBuilderName,muonBuilder);


      //Adams additions for KINK
      std::vector<Trajectory> refitted=theGlbRefitter->refit(*refL3,1);
      if(refitted.size()>0) {
	Trajectory muon = refitted.front();
	typedef TransientTrackingRecHit::ConstRecHitPointer 	ConstRecHitPointer;
	typedef ConstRecHitPointer RecHit;
	typedef std::vector<TrajectoryMeasurement>::const_iterator TMI;
	
	std::vector<TrajectoryMeasurement> meas = muon.measurements();

	for ( TMI m = meas.begin(); m != meas.end(); m++ ) {
	  TransientTrackingRecHit::ConstRecHitPointer hit = m->recHit();
	  RecHit rhit = (*m).recHit();

	  double phi2 = rhit->globalPosition().phi();
	  if ( phi2 < 0 ) phi2 = 2*M_PI + phi2;

	  GlobalPoint hitPos = rhit->globalPosition();	  
	  GlobalError hitErr = rhit->globalPositionError();
	  double error = hitErr.phierr(hitPos);  // error squared
	  
	  (*xForThisL3TM).push_back(rhit->globalPosition().x());
	  (*yForThisL3TM).push_back(rhit->globalPosition().y());
	  (*zForThisL3TM).push_back(rhit->globalPosition().z());
	  (*phiForThisL3TM).push_back(phi2);
	  (*errorForThisL3TM).push_back(error);

	  const TrajectoryStateOnSurface& tsos = (*m).predictedState();
	  if ( tsos.isValid() && rhit->isValid() && rhit->hit()->isValid()
	       && !std::isinf(rhit->localPositionError().xx()) //this is paranoia induced by reported case
	       && !std::isinf(rhit->localPositionError().xy()) //it's better to track down the origin of bad numbers
	       && !std::isinf(rhit->localPositionError().yy())
	       //	       && !std::isnan(rhit->localPositionError().xx()) //this is paranoia induced by reported case
	       //	       && !std::isnan(rhit->localPositionError().xy()) //it's better to track down the origin of bad numbers
	       //	       && !std::isnan(rhit->localPositionError().yy())
	       ) {
	    
	    double phi1 = tsos.globalPosition().phi();
	    if ( phi1 < 0 ) phi1 = 2*M_PI + phi1;
	    
	    (*xForThisL3TSOS).push_back(tsos.globalPosition().x());
	    (*yForThisL3TSOS).push_back(tsos.globalPosition().y());
	    (*zForThisL3TSOS).push_back(tsos.globalPosition().z());
	    (*phiForThisL3TSOS).push_back(phi1);

	    double phi2 = rhit->globalPosition().phi();
	    if ( phi2 < 0 ) phi2 = 2*M_PI + phi2;
	    
	    double diff = fabs(phi1 - phi2);
	    if ( diff > M_PI ) diff = 2*M_PI - diff;
	    
	    GlobalPoint hitPos = rhit->globalPosition();
	    
	    GlobalError hitErr = rhit->globalPositionError();
	    //LogDebug(theCategory)<<"hitPos " << hitPos;
	    double error = hitErr.phierr(hitPos);  // error squared
	    double s = ( error > 0.0 ) ? (diff*diff)/error : (diff*diff);
	  }
	}	
      }
      //end Adams KINK 

      for (trackingRecHit_iterator l3Hit = refL3->recHitsBegin(); l3Hit != refL3->recHitsEnd(); ++l3Hit) {
	if ((*l3Hit)->isValid()) {
	  (*idsForThisL3).push_back((*l3Hit)->geographicalId().rawId());
	  (*subidsForThisL3).push_back((*l3Hit)->geographicalId().subdetId());
	  (*detsForThisL3).push_back((*l3Hit)->geographicalId().det());
	  (*statusForThisL3).push_back((*l3Hit)->type());
	  if ((*l3Hit)->geographicalId().det() == 1) { // Tracker
	    TransientTrackingRecHit::RecHitPointer globL3 = trackBuilder->build(&**l3Hit);
	    (*xForThisL3).push_back(globL3->globalPosition().x());
	    (*yForThisL3).push_back(globL3->globalPosition().y());
	    (*zForThisL3).push_back(globL3->globalPosition().z());
	  }
	  else if ((*l3Hit)->geographicalId().det() == 2) { // Muon System
	    nMuHitsForThisL3++;
	    TransientTrackingRecHit::RecHitPointer globL3 = muonBuilder->build(&**l3Hit);
	    (*xForThisL3).push_back(globL3->globalPosition().x());
	    (*yForThisL3).push_back(globL3->globalPosition().y());
	    (*zForThisL3).push_back(globL3->globalPosition().z());
	    // number station goes here
	    if ( (*l3Hit)->geographicalId().subdetId() == 1) { // DT hit 
	      const DTChamberId& id = DTChamberId((*l3Hit)->geographicalId());
	      (*stationsForThisL3).push_back(id.station());
	    }
	    if ( (*l3Hit)->geographicalId().subdetId() == 2) { // CSC hit
	      const CSCDetId& id = CSCDetId((*l3Hit)->geographicalId());
	      (*stationsForThisL3).push_back(id.station());
	    }
	    if ( (*l3Hit)->geographicalId().subdetId() == 3) { // RPC hit
	      const RPCDetId& id = RPCDetId((*l3Hit)->geographicalId());
	      (*stationsForThisL3).push_back(id.station());
	    }
	  }
	}
      }

      (*l3DetIds).insert(std::make_pair(iMu,*idsForThisL3));    
      (*l3SubdetIds).insert(std::make_pair(iMu,*subidsForThisL3));
      (*l3Component).insert(std::make_pair(iMu,*detsForThisL3));
      (*l3RecHitsStatus).insert(std::make_pair(iMu,*statusForThisL3));
      (*l3NMuHits).insert(std::make_pair(iMu,nMuHitsForThisL3));
      (*l3MuStationNumber).insert(std::make_pair(iMu,*stationsForThisL3));
      (*l3RecHitsX).insert(std::make_pair(iMu,*xForThisL3));
      (*l3RecHitsY).insert(std::make_pair(iMu,*yForThisL3));
      (*l3RecHitsZ).insert(std::make_pair(iMu,*zForThisL3));
      (*l3RecHitsXTM).insert(std::make_pair(iMu,*xForThisL3TM));
      (*l3RecHitsYTM).insert(std::make_pair(iMu,*yForThisL3TM));
      (*l3RecHitsZTM).insert(std::make_pair(iMu,*zForThisL3TM));
      (*l3RecHitsXTSOS).insert(std::make_pair(iMu,*xForThisL3TSOS));
      (*l3RecHitsYTSOS).insert(std::make_pair(iMu,*yForThisL3TSOS));
      (*l3RecHitsZTSOS).insert(std::make_pair(iMu,*zForThisL3TSOS));
      (*l3RecHitsErrorTM).insert(std::make_pair(iMu,*errorForThisL3TM));
      (*l3RecHitsPhiTM).insert(std::make_pair(iMu,*phiForThisL3TM));
      (*l3RecHitsPhiTSOS).insert(std::make_pair(iMu,*phiForThisL3TSOS));

      idsForThisL3->clear();
      subidsForThisL3->clear();
      detsForThisL3->clear();
      statusForThisL3->clear();
      stationsForThisL3->clear();
      xForThisL3->clear();
      yForThisL3->clear();
      zForThisL3->clear();
      xForThisL3TM->clear();
      yForThisL3TM->clear();
      zForThisL3TM->clear();
      xForThisL3TSOS->clear();
      yForThisL3TSOS->clear();
      zForThisL3TSOS->clear();
      errorForThisL3TM->clear();
      phiForThisL3TM->clear();
      phiForThisL3TSOS->clear();
      nMuHitsForThisL3 = 0;

      /*  
      // Find the correct L2 muon using track links.  Start with borrowed code.
      bool correctTrackLink = false;
      uint matchingLink=0;
      for(uint iLink=0;iLink!=l3ToL2Links->size();iLink++){
      // L3 is treated as a global track.
      if ((*l3ToL2Links)[iLink].globalTrack() == refL3){
      //this L3muon is coming from this L2 track
      correctTrackLink =true;
      matchingLink=iLink;
      }//same track ref
      }//loop over track links
      if(correctTrackLink) {
      reco::TrackRef refL2FromTrackLinks = (*l3ToL2Links)[matchingLink].standAloneTrack();
      // Get match indices here
      for(int iL2 = 0;iL2!=nL2;++iL2){
      reco::TrackRef refL2ToMatch(l2Muons, iL2);
      if (refL2FromTrackLinks == refL2ToMatch) {
      // Fill these so we know which L2 muon seeded this L3 muon
      (*indexL2SeedingL3).push_back(iL2);
      (*indexL3SeededFromL2).push_back(iL3);
      }
      }
      }
      */
      
      // get the CAL deposits associated with this muon
      reco::IsoDeposit calDeposit = caloDepositExtractor->deposit(iEvent, iSetup, *refL3);
      // cutting for the L3 muon cal isolation (just getting the cuts and vetos)
      muonisolation::Cuts::CutSpec calo_cuts_here = L2IsoCalCuts(refL3->eta());
      // and deposit for the L3 muon cal isolation
      double conesize = calo_cuts_here.conesize;
      (*l3CalIsoDeposit).push_back(calDeposit.depositWithin(conesize));
      // now for the L3 tracking things
      reco::IsoDeposit trackDeposit = trackDepositExtractor->deposit(iEvent, iSetup, *refL3);
      reco::IsoDeposit::const_iterator pixEnd = trackDeposit.end();
      double dRMinDelPt = 1;
      double dRMin = 1;
      double DelPtMin = 999;
      for (reco::IsoDeposit::const_iterator pix = trackDeposit.begin(); pix != pixEnd; ++pix) {
	double tempDelphi = pix.phi() - refL3->phi();
	if (tempDelphi > TMath::Pi()) tempDelphi = (2 * TMath::Pi()) - tempDelphi;
	double tempdR = sqrt(((pix.eta() - refL3->eta()) * (pix.eta() - refL3->eta())) + (tempDelphi * tempDelphi));
	double tempPt = pix.value();
	if (tempdR < dRMin) dRMin = tempdR;
	if (fabs(refL3->pt() - tempPt) < DelPtMin) {
	  dRMinDelPt = tempdR;
	  DelPtMin = fabs(refL3->pt() - tempPt);
	}
      }
      (*l3IsoTrackDR).push_back(dRMin);
      (*l3IsoTrackDRMinDelPt).push_back(dRMinDelPt);
      trackVetos.push_back(trackDeposit.veto());
      const muonisolation::Cuts::CutSpec & trackCut = L3IsoTrackCuts(refL3->eta());
      // get the Tracking deposits for our muon
      (*l3TrackIsoDeposit).push_back(trackDeposit.depositWithin(trackCut.conesize, trackVetos));

      bool associated = false;
      // With the detector-level things filled, time to start doing the associations to sim
      int sim_index = 0;
      TrackingParticleRef muTP;
      TrackingParticleRef anyTP;
      for (reco::RecoToSimCollection::const_iterator findRefL3 = l3RecSimColl.begin(); findRefL3 != l3RecSimColl.end(); ++findRefL3) {
	const edm::RefToBase<reco::Track> & l3RecSimMatch = findRefL3->key;
	if (l3RecSimMatch->pt() == refL3->pt()) {
	  associated = true;
	  const std::vector<std::pair<TrackingParticleRef,double> > & tp = findRefL3->val;
	  const TrackingParticleRef & trp = tp.begin()->first;
	  
	  //Start adam hack to only use muons

	  //end adam hack to only use muons

	  (*l3AssociationVar).push_back(tp.begin()->second);
	  
	  int particle_ID = trp->pdgId();
	  //	int myBin = wantMotherBin.GetBinNum(particle_ID);
	  (*l3AssociationPdgId).push_back(particle_ID);
	  LogDebug("SpecialBit")<<"SpecialBit pdgId " << particle_ID;
	  unsigned int thisBit = getBit(trp);
	  //
	  int muClass = 0;
	  if (isPrimaryMuon(thisBit)) muClass = 1;
	  else if (isSiliconMuon(thisBit)) muClass = 2;
	  else if (isCalConversionMuon(thisBit)) muClass = 3;
	  else if (isOtherMuon(thisBit)) muClass = 4;
	  else muClass = 5;
	  LogDebug("SpecialBit")<<"RecoTree muon " << iMu << " of " << muonColl.size() << " has theBit: " << thisBit << " and muClass " << muClass << " and pT " << refL3->pt();
	  //
	  (*l3AssociationMyBit).push_back(thisBit);
	  (*l3AssociationVtxX).push_back(trp->vertex().x());
	  (*l3AssociationVtxY).push_back(trp->vertex().y());
	  (*l3AssociationVtxZ).push_back(trp->vertex().z());
	  
	  if(abs(particle_ID) == 13){
	    //Adam 1 start
	    int iSimMu = -1;
	    for (unsigned int iSim = 0; iSim != (*TPtracks).size(); iSim++) {
	      TrackingParticleRef trp2(TPtracks, iSim);
	      int particle_ID2 = trp2->pdgId();
	      if(abs(particle_ID2) != 13) continue;
	      iSimMu++;
	      (*l3AssociatedSimMuonIndex).push_back(iSimMu);
	    }
	    //Adam 1 end
	    // put in the associated pt,eta,phi
	    (*l3AssociatedSimMuonPt).push_back(trp->pt());
	    (*l3AssociatedSimMuonEta).push_back(trp->eta());
	    (*l3AssociatedSimMuonPhi).push_back(trp->phi());
	    if (fabs(trp->phi() - refL3->phi()) > 1) {
	      //Note: keeping this in.  This happens sometimes when the associator used is the 
	      //steppingHelixPropagatorAny
	      edm::LogInfo("MuonRecoTreeUtility") << "Something's gone wrong here. First our indexes";
	      edm::LogInfo("MuonRecoTreeUtility") << "iL3, sim_index = " << iMu <<" " << sim_index;
	      edm::LogInfo("MuonRecoTreeUtility") << "What about recSimMatch vs trp phi?" << l3RecSimMatch->phi() <<" " << trp->phi();
	    }
	    // put in the detIDs for this sim muon
	    std::vector<int> *idsForSimL3 = new std::vector<int>;
	    std::vector<int> *stationsForSim = new std::vector<int>;
	    int simHitCounter = 0;
	    int simMuHitCounter = 0;
	    for (PSimHitContainer::const_iterator l3SimHit = trp->pSimHit_begin(); l3SimHit != trp->pSimHit_end(); ++l3SimHit) {
	      (*idsForSimL3).push_back((*l3SimHit).detUnitId());
	      DetId theDetUnitId(l3SimHit->detUnitId());
	      int detector = theDetUnitId.det();
	      int subdetector = theDetUnitId.subdetId();
	      if (detector == 2) { //Muon system
		simMuHitCounter ++;
		if (subdetector == 1) { //DT
		  const DTChamberId& id = DTChamberId(l3SimHit->detUnitId());
		  (*stationsForSim).push_back(id.station());
		}
		if (subdetector == 2) { //CSC
		  const CSCDetId& id=CSCDetId(l3SimHit->detUnitId());
		  (*stationsForSim).push_back(id.station());
		}
		if (subdetector == 3) { //RPC
		  const RPCDetId& id = RPCDetId(l3SimHit->detUnitId());
		  (*stationsForSim).push_back(id.station());
		}
	      }
	      simHitCounter++;
	    }
	    (*l3AssociatedSimMuonNHits).push_back(simHitCounter);
	    (*l3AssociatedSimMuonDetIds).insert(std::make_pair(sim_index,*idsForSimL3));
	    (*l3AssociatedSimMuonMuStationNumber).insert(std::make_pair(sim_index,*stationsForSim));
	    (*l3AssociatedSimMuonNMuHits).insert(std::make_pair(sim_index,simMuHitCounter));
	    idsForSimL3->clear();
	    stationsForSim->clear();
	    sim_index++;
	    //---------------------- MOTHERHOOD ----------------------------
	    //find the parent of tracking particle
	    for(TrackingParticle::g4t_iterator isimtk = trp->g4Track_begin();isimtk!=trp->g4Track_end();isimtk++)
	      {
		LogDebug(theCategory)<<"I am here 1";
		if(isimtk->type()==13||isimtk->type()==-13)
		  {
		    // This is the sim track for this tracking particle.  Time to put in the parameters
		    FreeTrajectoryState 
		      ftsAtProduction(GlobalPoint(trp->vertex().x(),trp->vertex().y(),trp->vertex().z()),
				      GlobalVector(isimtk->momentum().x(),isimtk->momentum().y(),isimtk->momentum().z()),
				      TrackCharge(trp->charge()),
				      field.product());
		    TSCBLBuilderNoMaterial tscblBuilder;
		    TrajectoryStateClosestToBeamLine tsAtClosestApproach = tscblBuilder(ftsAtProduction,bs);//as in TrackProducerAlgorithm
		    GlobalPoint v1 = tsAtClosestApproach.isValid() ? tsAtClosestApproach.trackStateAtPCA().position() : GlobalPoint() ;
		    GlobalVector p = tsAtClosestApproach.isValid() ? tsAtClosestApproach.trackStateAtPCA().momentum() : GlobalVector() ;
		    GlobalPoint v(v1.x()-bs.x0(),v1.y()-bs.y0(),v1.z()-bs.z0());
		    
		    double qoverpSim = tsAtClosestApproach.isValid() ? tsAtClosestApproach.trackStateAtPCA().charge()/p.mag() : 0.;
		    double lambdaSim = M_PI/2-p.theta();
		    double dxySim    = (-v.x()*sin(p.phi())+v.y()*cos(p.phi()));
		    double dzSim     = v.z() - (v.x()*p.x()+v.y()*p.y())/p.perp() * p.z()/p.perp();
		    
		    (*l3AssociatedSimMuonDsz).push_back(dzSim);
		    (*l3AssociatedSimMuonDxy).push_back(dxySim);
		    (*l3AssociatedSimMuonLambda).push_back(lambdaSim);
		    (*l3AssociatedSimMuonQoverP).push_back(qoverpSim);
		    
		    //calculate mother hood
		    MotherSearch mother(&*isimtk, SimTk, SimVtx, hepmc);
		    //FIXME, use reco::Particle mother.mother();
		    //                double pt,eta,phi;
		    //                int parentID;
		    //                int motherBinNumber;
		    
		    if (mother.IsValid()){
		      if (mother.SimIsValid()){
			(*l3ParentID).push_back(mother.Sim_mother->type());
			(*l3MotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Sim_mother->type()));
		      }
		      else {
			(*l3ParentID).push_back(mother.Gen_mother->pdg_id());
			(*l3MotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id()));
		      }
		      //do it once per tracking particle once it succeeds
		      break;
		    }
		    else{
		      // This handles cases when we have an associated sim muon, but "tricky" muon without 
		      // valid parent (e.g. singleMu)
		      (*l3ParentID).push_back(0);
		      (*l3MotherBinNumber).push_back(wantMotherBin.GetBinNum(0));
		      edm::LogError(theCategory)<<"tricky muon from TrackingParticle.";
		    }
		  }//sim track is a muon
		else{
		  (*l3ParentID).push_back(isimtk->type());
		  (*l3MotherBinNumber).push_back(777);
		  edm::LogError(theCategory)<<"the sim track attached to the tracking particle is not a muon.";
		}
	      }  //loop over SimTrack of tracking particle
	  }// particle_ID == 13
	  else{
	    //a reco muon is associated to something else than a muon
	    edm::LogError(theCategory)<<"a reconstructed muon is associated to: "<<particle_ID;
	    (*l3ParentID).push_back(particle_ID);
	    (*l3MotherBinNumber).push_back(-777);
	  }
	}//track has an association
	else{
	  //this track was not associated.
	  edm::LogError(theCategory)<<"a reconstructed muon is not associated to a muonTP";
	}
      }
      if (associated) {
	//      std::cout << "Associated..." << std::endl;
	(*l3IsAssociated).push_back(1);
      }
      else {
	(*l3IsAssociated).push_back(0);
	(*l3AssociationVar).push_back(-999);
	(*l3AssociatedSimMuonIndex).push_back(-999);
	(*l3AssociatedSimMuonPt).push_back(-999);
	(*l3AssociatedSimMuonEta).push_back(-999);
	(*l3AssociatedSimMuonPhi).push_back(-999);
	(*l3AssociatedSimMuonNHits).push_back(-999);
	(*l3AssociatedSimMuonQoverP).push_back(-999);
	(*l3AssociatedSimMuonLambda).push_back(-999);
	(*l3AssociatedSimMuonDxy).push_back(-999);
	(*l3AssociatedSimMuonDsz).push_back(-999);
	(*l3ParentID).push_back(-666);
	(*l3MotherBinNumber).push_back(-999);
      }
    } else { //ifGlbTrack isAvailable  //loop over l3Muons
      double crap = -999;
      (*l3P).push_back(crap);
      (*l3Px).push_back(crap);
      (*l3Py).push_back(crap);
      (*l3Pz).push_back(crap);
      (*l3Pt).push_back(crap);
      (*l3PtError).push_back(crap);
      (*l3Pt90).push_back(crap);
      (*l3Eta).push_back(crap);
      (*l3EtaError).push_back(crap);
      (*l3Phi).push_back(crap);
      (*l3PhiError).push_back(crap);
      (*l3D0).push_back(crap);
      (*l3D0Error).push_back(crap);
      (*l3NHits).push_back(0);
      (*l3Charge).push_back(crap);
      (*l3Chi2).push_back(crap);
      (*l3Ndof).push_back(crap);
      // Fill in the track fitting parameters (with phi filled in above)
      (*l3Dsz).push_back(crap);
      (*l3DszError).push_back(crap);
      (*l3Dxy).push_back(crap);
      (*l3DxyError).push_back(crap);
      (*l3Lambda).push_back(crap);
      (*l3LambdaError).push_back(crap);
      (*l3Qoverp).push_back(crap);
      (*l3QoverpError).push_back(crap);
      //adam (*l3ErrorMatrix).push_back(NULL);
      
      (*l3DetIds).insert(std::make_pair(iMu,0));    
      (*l3SubdetIds).insert(std::make_pair(iMu,0));
      (*l3Component).insert(std::make_pair(iMu,0));
      (*l3RecHitsStatus).insert(std::make_pair(iMu,0));
      (*l3NMuHits).insert(std::make_pair(iMu,0));
      (*l3MuStationNumber).insert(std::make_pair(iMu,0));
      (*l3RecHitsX).insert(std::make_pair(iMu,0));
      (*l3RecHitsY).insert(std::make_pair(iMu,0));
      (*l3RecHitsZ).insert(std::make_pair(iMu,0));
      (*l3RecHitsXTM).insert(std::make_pair(iMu,0));
      (*l3RecHitsYTM).insert(std::make_pair(iMu,0));
      (*l3RecHitsZTM).insert(std::make_pair(iMu,0));
      (*l3RecHitsXTSOS).insert(std::make_pair(iMu,0));
      (*l3RecHitsYTSOS).insert(std::make_pair(iMu,0));
      (*l3RecHitsZTSOS).insert(std::make_pair(iMu,0));
      (*l3RecHitsPhiTM).insert(std::make_pair(iMu,0));
      (*l3RecHitsErrorTM).insert(std::make_pair(iMu,0));
      (*l3RecHitsPhiTSOS).insert(std::make_pair(iMu,0));
      
      (*l3IsoTrackDR).push_back(crap);
      (*l3IsoTrackDRMinDelPt).push_back(crap);
      (*l3TrackIsoDeposit).push_back(0);

      (*l3IsAssociated).push_back(0);      
      (*l3AssociationVar).push_back(crap);
      (*l3AssociationPdgId).push_back(-999);
      (*l3AssociationMyBit).push_back(0);
      (*l3AssociationVtxX).push_back(crap);
      (*l3AssociationVtxY).push_back(crap);
      (*l3AssociationVtxZ).push_back(crap);
      (*l3AssociatedSimMuonIndex).push_back(999);
      (*l3AssociatedSimMuonPt).push_back(crap);
      (*l3AssociatedSimMuonEta).push_back(crap);
      (*l3AssociatedSimMuonPhi).push_back(crap);
      
      (*l3AssociatedSimMuonNHits).push_back(0);
      (*l3AssociatedSimMuonDetIds).insert(std::make_pair(iMu,0));
      (*l3AssociatedSimMuonMuStationNumber).insert(std::make_pair(iMu,0));
      (*l3AssociatedSimMuonNMuHits).insert(std::make_pair(iMu,0));
      

      //      (*l3AssociationVar).push_back(-999);
      //      (*l3AssociatedSimMuonPt).push_back(-999);
      //      (*l3AssociatedSimMuonEta).push_back(-999);
      //      (*l3AssociatedSimMuonPhi).push_back(-999);
      //      (*l3AssociatedSimMuonNHits).push_back(-999);
      (*l3AssociatedSimMuonQoverP).push_back(-999);
      (*l3AssociatedSimMuonLambda).push_back(-999);
      (*l3AssociatedSimMuonDxy).push_back(-999);
      (*l3AssociatedSimMuonDsz).push_back(-999);
      (*l3ParentID).push_back(-666);
      (*l3MotherBinNumber).push_back(-999);
    }
    
    if(tkTrack.isAvailable()) {
      int iTk = iMu;  
      const reco::TrackRef refTk(tkTrack);
      
      // Fill in basic information of tkMuons
      (*tkTrackP).push_back(refTk->p());
      (*tkTrackPx).push_back(refTk->px());
      (*tkTrackPy).push_back(refTk->py());
      (*tkTrackPz).push_back(refTk->pz());
      (*tkTrackPt).push_back(refTk->pt());
      (*tkTrackPtError).push_back(refTk->ptError());
      (*tkTrackEta).push_back(refTk->eta());
      (*tkTrackEtaError).push_back(refTk->etaError());
      (*tkTrackPhi).push_back(refTk->phi());
      (*tkTrackPhiError).push_back(refTk->phiError());
      (*tkTrackD0).push_back(refTk->d0());
      (*tkTrackD0Error).push_back(refTk->d0Error());
      (*tkTrackNHits).push_back(refTk->recHitsSize());
      (*tkTrackCharge).push_back(refTk->charge());
      (*tkTrackChi2).push_back(refTk->chi2());
      (*tkTrackNdof).push_back(refTk->ndof());
      // Fill in the track fitting parameters (with phi filled in above)
      (*tkTrackDsz).push_back(refTk->dsz());
      (*tkTrackDszError).push_back(refTk->dszError());
      (*tkTrackDxy).push_back(refTk->dxy());
      (*tkTrackDxyError).push_back(refTk->dxyError());
      (*tkTrackLambda).push_back(refTk->lambda());
      (*tkTrackLambdaError).push_back(refTk->lambdaError());
      (*tkTrackQoverp).push_back(refTk->qoverp());
      (*tkTrackQoverpError).push_back(refTk->qoverpError());
      (*tkTrackErrorMatrix).push_back(refTk->covariance());
      
      std::vector<int> *idsForThisTk = new std::vector<int>;
      std::vector<int> *subidsForThisTk = new std::vector<int>;
      std::vector<int> *statusForThisTk = new std::vector<int>;
      std::vector<double> *xForThisTk = new std::vector<double>;
      std::vector<double> *yForThisTk = new std::vector<double>;
      std::vector<double> *zForThisTk = new std::vector<double>;
      
      edm::ESHandle<TransientTrackingRecHitBuilder> trackBuilder;
      std::string trackBuilderName = "WithTrackAngle";
      iSetup.get<TransientRecHitRecord>().get(trackBuilderName,trackBuilder);
      
      edm::LogInfo("MuonRecoTreeUtility") << "iterating over tracking rechits";
      for (trackingRecHit_iterator tkHit = refTk->recHitsBegin(); tkHit != refTk->recHitsEnd(); ++tkHit) {     
	if ((*tkHit)->isValid()) {
	  (*idsForThisTk).push_back((*tkHit)->geographicalId().rawId());
	  (*subidsForThisTk).push_back((*tkHit)->geographicalId().subdetId());
	  (*statusForThisTk).push_back((*tkHit)->type());
	  if ((*tkHit)->geographicalId().det() == 1) {
	    TransientTrackingRecHit::RecHitPointer globTk = trackBuilder->build(&**tkHit);
	    (*xForThisTk).push_back(globTk->globalPosition().x());
	    (*yForThisTk).push_back(globTk->globalPosition().y());
	    (*zForThisTk).push_back(globTk->globalPosition().z());
	  }
	  else edm::LogInfo("MuonRecoTreeUtility") << "rechit found in detector " << (*tkHit)->geographicalId().det();
	}
      }
      edm::LogInfo("MuonRecoTreeUtility") << "iteration finished";
      
      (*tkTrackDetIds).insert(std::make_pair(iMu,*idsForThisTk));
      (*tkTrackSubdetIds).insert(std::make_pair(iMu,*subidsForThisTk));
      (*tkTrackRecHitsStatus).insert(std::make_pair(iMu,*statusForThisTk));
      (*tkTrackRecHitsX).insert(std::make_pair(iMu,*xForThisTk));
      (*tkTrackRecHitsY).insert(std::make_pair(iMu,*yForThisTk));
      (*tkTrackRecHitsZ).insert(std::make_pair(iMu,*zForThisTk));
      
      idsForThisTk->clear();
      subidsForThisTk->clear();
      statusForThisTk->clear();
      xForThisTk->clear();
      yForThisTk->clear();
      zForThisTk->clear();
      
      bool associated = false;
      // With the detector-level things filled, time to start doing the associations to sim
      int sim_index = 0;
      for (reco::RecoToSimCollection::const_iterator findRefTk = tkRecSimColl.begin(); findRefTk != tkRecSimColl.end(); ++findRefTk) {
	const edm::RefToBase<reco::Track> & tkRecSimMatch = findRefTk->key;
	if (tkRecSimMatch->pt() == refTk->pt()) {
	  associated = true;
	  const std::vector<std::pair<TrackingParticleRef,double> > & tp = findRefTk->val;
	  const TrackingParticleRef & trp = tp.begin()->first;
	  
	  (*tkTrackAssociationVar).push_back(tp.begin()->second);
	  
	  int particle_ID = trp->pdgId();
	  //	int myBin = wantMotherBin.GetBinNum(particle_ID);
	  
	  if(abs(particle_ID) == 13){
	    // put in the associated pt,eta,phi
	    (*tkTrackAssociatedSimMuonPt).push_back(trp->pt());
	    (*tkTrackAssociatedSimMuonEta).push_back(trp->eta());
	    (*tkTrackAssociatedSimMuonPhi).push_back(trp->phi());
	    if (fabs(trp->phi() - refTk->phi()) > 1) {
	      //Note: keeping this in.  This happens sometimes when the associator used is the
	      //steppingHelixPropagatorAny
	      edm::LogInfo("MuonRecoTreeUtility") << "Something's gone wrong here. First our indexes";
	      edm::LogInfo("MuonRecoTreeUtility") << "iTk, sim_index = " << iMu <<" " << sim_index;
	      edm::LogInfo("MuonRecoTreeUtility") << "What about recSimMatch vs trp phi?" << tkRecSimMatch->phi() <<" " << trp->phi();
	    }
	    // put in the detIDs for this sim muon
	    std::vector<int> *idsForSimTk = new std::vector<int>;
	    int simHitCounter = 0;
	    for (PSimHitContainer::const_iterator tkSimHit = trp->pSimHit_begin(); tkSimHit != trp->pSimHit_end(); ++tkSimHit) {
	      (*idsForSimTk).push_back((*tkSimHit).detUnitId());
	      DetId theDetUnitId(tkSimHit->detUnitId());
	      simHitCounter++;
	    }
	    (*tkTrackAssociatedSimMuonNHits).push_back(simHitCounter);
	    (*tkTrackAssociatedSimMuonDetIds).insert(std::make_pair(sim_index,*idsForSimTk));
	    idsForSimTk->clear();
	    sim_index++;
	    //---------------------- MOTHERHOOD ---------------------------
	    //find the parent of tracking particle
	    for(TrackingParticle::g4t_iterator isimtk = trp->g4Track_begin();isimtk!=trp->g4Track_end();isimtk++)
	      {
		LogDebug(theCategory)<<"I am here 1";
		if(isimtk->type()==13||isimtk->type()==-13)
		  {
		    // This is the sim track for this tracking particle.  Time to put in the parameters
		    FreeTrajectoryState
		      ftsAtProduction(GlobalPoint(trp->vertex().x(),trp->vertex().y(),trp->vertex().z()),
				      GlobalVector(isimtk->momentum().x(),isimtk->momentum().y(),isimtk->momentum().z()),
				      TrackCharge(trp->charge()),
				      field.product());
		    TSCBLBuilderNoMaterial tscblBuilder;
		    TrajectoryStateClosestToBeamLine tsAtClosestApproach = tscblBuilder(ftsAtProduction,bs);//as in TrackProducerAlgorithm
		    GlobalPoint v1 = tsAtClosestApproach.isValid() ? tsAtClosestApproach.trackStateAtPCA().position() : GlobalPoint();
		    GlobalVector p = tsAtClosestApproach.isValid() ? tsAtClosestApproach.trackStateAtPCA().momentum() : GlobalVector();
		    GlobalPoint v(v1.x()-bs.x0(),v1.y()-bs.y0(),v1.z()-bs.z0());
		    
		    double qoverpSim = tsAtClosestApproach.isValid() ? tsAtClosestApproach.trackStateAtPCA().charge()/p.mag() : 0.;
		    double lambdaSim = M_PI/2-p.theta();
		    double dxySim    = (-v.x()*sin(p.phi())+v.y()*cos(p.phi()));
		    double dzSim     = v.z() - (v.x()*p.x()+v.y()*p.y())/p.perp() * p.z()/p.perp();
		    
		    (*tkTrackAssociatedSimMuonDsz).push_back(dzSim);
		    (*tkTrackAssociatedSimMuonDxy).push_back(dxySim);
		    (*tkTrackAssociatedSimMuonLambda).push_back(lambdaSim);
		    (*tkTrackAssociatedSimMuonQoverP).push_back(qoverpSim);
		    
		    //calculate mother hood
		    MotherSearch mother(&*isimtk, SimTk, SimVtx, hepmc);
		    //FIXME, use reco::Particle mother.mother();
		    //                double pt,eta,phi;
		    //                int parentID;
		    //                int motherBinNumber;
		    
		    if (mother.IsValid()){
		      if (mother.SimIsValid()){
			(*tkTrackParentID).push_back(mother.Sim_mother->type());
			(*tkTrackMotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Sim_mother->type()));
		      }
		      else {
			(*tkTrackParentID).push_back(mother.Gen_mother->pdg_id());
			(*tkTrackMotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id()));
		      }
		      //do it once per tracking particle once it succeeds
		      break;
		    }
		    else{
		      // This handles cases when we have an associated sim muon, but "tricky" muon without
		      // valid parent (e.g. singleMu)
		      (*tkTrackParentID).push_back(0);
		      (*tkTrackMotherBinNumber).push_back(wantMotherBin.GetBinNum(0));
		      edm::LogError(theCategory)<<"tricky muon from TrackingParticle.";
		    }
		  }//sim track is a muon
		else{
		  (*tkTrackParentID).push_back(isimtk->type());
		  (*tkTrackMotherBinNumber).push_back(777);
		  edm::LogError(theCategory)<<"the sim track attached to the tracking particle is not a muon.";
		}
	      }//loop over SimTrack of tracking particle
	  }//muon associated
	  else{
	    //a reco muon is associated to something else than a muon
	    edm::LogError(theCategory)<<"a reconstructed muon is associated to: "<<particle_ID;
	    (*tkTrackParentID).push_back(particle_ID);
	    (*tkTrackMotherBinNumber).push_back(-777);
	  }
	}//track has an association
	else{
	  //this track was not associated.
	  edm::LogError(theCategory)<<"a reconstructed muon is not associated.";
	}
      }
      if (associated) {
	//      std::cout << "Associated..." << std::endl;
	(*tkTrackIsAssociated).push_back(1);
      }
      else {
	(*tkTrackIsAssociated).push_back(0);
	(*tkTrackAssociationVar).push_back(-999);
	(*tkTrackAssociatedSimMuonPt).push_back(-999);
	(*tkTrackAssociatedSimMuonEta).push_back(-999);
	(*tkTrackAssociatedSimMuonPhi).push_back(-999);
	(*tkTrackAssociatedSimMuonNHits).push_back(-999);
	(*tkTrackAssociatedSimMuonQoverP).push_back(-999);
	(*tkTrackAssociatedSimMuonLambda).push_back(-999);
	(*tkTrackAssociatedSimMuonDxy).push_back(-999);
	(*tkTrackAssociatedSimMuonDsz).push_back(-999);
	(*tkTrackParentID).push_back(-666);
	(*tkTrackMotherBinNumber).push_back(-999);
      }
    } else { // loop over hltL3TkTracksFromL2
      double crap = -999;
      (*tkTrackP).push_back(crap);
      (*tkTrackPx).push_back(crap);
      (*tkTrackPy).push_back(crap);
      (*tkTrackPz).push_back(crap);
      (*tkTrackPt).push_back(crap);
      (*tkTrackPtError).push_back(crap);
      //(*tkTrackPt90).push_back(crap);
      (*tkTrackEta).push_back(crap);
      (*tkTrackEtaError).push_back(crap);
      (*tkTrackPhi).push_back(crap);
      (*tkTrackPhiError).push_back(crap);
      (*tkTrackD0).push_back(crap);
      (*tkTrackD0Error).push_back(crap);
      (*tkTrackNHits).push_back(0);
      (*tkTrackCharge).push_back(crap);
      (*tkTrackChi2).push_back(crap);
      (*tkTrackNdof).push_back(crap);
      // Fill in the track fitting parameters (with phi filled in above)
      (*tkTrackDsz).push_back(crap);
      (*tkTrackDszError).push_back(crap);
      (*tkTrackDxy).push_back(crap);
      (*tkTrackDxyError).push_back(crap);
      (*tkTrackLambda).push_back(crap);
      (*tkTrackLambdaError).push_back(crap);
      (*tkTrackQoverp).push_back(crap);
      (*tkTrackQoverpError).push_back(crap);
      //adam (*l3ErrorMatrix).push_back(NULL);
      
      (*tkTrackDetIds).insert(std::make_pair(iMu,0));    
      (*tkTrackSubdetIds).insert(std::make_pair(iMu,0));
      //(*tkTrackComponent).insert(std::make_pair(iMu,0));
      (*tkTrackRecHitsStatus).insert(std::make_pair(iMu,0));
      //(*tkTrackNMuHits).insert(std::make_pair(iMu,0));
      //(*tkTrackMuStationNumber).insert(std::make_pair(iMu,0));
      (*tkTrackRecHitsX).insert(std::make_pair(iMu,0));
      (*tkTrackRecHitsY).insert(std::make_pair(iMu,0));
      (*tkTrackRecHitsZ).insert(std::make_pair(iMu,0));
      
      //(*tkTrackIsoTrackDR).push_back(crap);
      //(*tkTrackIsoTrackDRMinDelPt).push_back(crap);
      //(*tkTrackTrackIsoDeposit).push_back(0);

      (*tkTrackIsAssociated).push_back(0);      
      (*tkTrackAssociationVar).push_back(crap);
      (*tkTrackAssociationPdgId).push_back(-999);
      (*tkTrackAssociationMyBit).push_back(0);
      (*tkTrackAssociationVtxX).push_back(crap);
      (*tkTrackAssociationVtxY).push_back(crap);
      (*tkTrackAssociationVtxZ).push_back(crap);
      (*tkTrackAssociatedSimMuonPt).push_back(crap);
      (*tkTrackAssociatedSimMuonEta).push_back(crap);
      (*tkTrackAssociatedSimMuonPhi).push_back(crap);
      
      (*tkTrackAssociatedSimMuonNHits).push_back(0);
      (*tkTrackAssociatedSimMuonDetIds).insert(std::make_pair(iMu,0));
      //(*tkTrackAssociatedSimMuonMuStationNumber).insert(std::make_pair(iMu,0));
      //(*tkTrackAssociatedSimMuonNMuHits).insert(std::make_pair(iMu,0));
      
      //      (*tkTrackAssociationVar).push_back(-999);
      //      (*tkTrackAssociatedSimMuonPt).push_back(-999);
      //      (*tkTrackAssociatedSimMuonEta).push_back(-999);
      //      (*tkTrackAssociatedSimMuonPhi).push_back(-999);
      //      (*tkTrackAssociatedSimMuonNHits).push_back(-999);
      (*tkTrackAssociatedSimMuonQoverP).push_back(-999);
      (*tkTrackAssociatedSimMuonLambda).push_back(-999);
      (*tkTrackAssociatedSimMuonDxy).push_back(-999);
      (*tkTrackAssociatedSimMuonDsz).push_back(-999);
      (*tkTrackParentID).push_back(-666);
      (*tkTrackMotherBinNumber).push_back(-999);
    }
    
    if(staTrack.isAvailable()) {  
      const reco::TrackRef refL2(staTrack);
      int iL2 = iMu;

      // Fill in basic information of l2Muons
      (*l2P).push_back(refL2->p());
      (*l2Px).push_back(refL2->px());
      (*l2Py).push_back(refL2->py());
      (*l2Pz).push_back(refL2->pz());
      (*l2Pt).push_back(refL2->pt());
      (*l2PtError).push_back(refL2->ptError());
      (*l2Pt90).push_back(pt90(refL2,iEvent));
      (*l2Eta).push_back(refL2->eta());
      (*l2EtaError).push_back(refL2->etaError());
      (*l2Phi).push_back(refL2->phi());
      (*l2PhiError).push_back(refL2->phiError());
      (*l2D0).push_back(refL2->d0());
      (*l2D0Error).push_back(refL2->d0Error());
      (*l2NHits).push_back(refL2->recHitsSize());
      (*l2Charge).push_back(refL2->charge());
      (*l2Chi2).push_back(refL2->chi2());
      (*l2Ndof).push_back(refL2->ndof());
      // Fill in the track fitting parameters (with phi filled in above)
      (*l2Dsz).push_back(refL2->dsz());
      (*l2DszError).push_back(refL2->dszError());
      (*l2Dxy).push_back(refL2->dxy());
      (*l2DxyError).push_back(refL2->dxyError());
      (*l2Lambda).push_back(refL2->lambda());
      (*l2LambdaError).push_back(refL2->lambdaError());
      (*l2Qoverp).push_back(refL2->qoverp());
      (*l2QoverpError).push_back(refL2->qoverpError());
      (*l2ErrorMatrix).push_back(refL2->covariance());
      // Filling in THE muon error matrix
      //    edm::LogInfo("MuonRecoTreeUtility") << "Trying to fill the muon error matrix.";
      
      std::vector<double> *matrixValuesForThisL2 = new std::vector<double>;
      std::vector<int> *idsForThisL2 = new std::vector<int>;
      std::vector<int> *subidsForThisL2 = new std::vector<int>;
      std::vector<int> *detsForThisL2 = new std::vector<int>;
      std::vector<int> *statusForThisL2 = new std::vector<int>;
      std::vector<double> *xForThisL2 = new std::vector<double>;
      std::vector<double> *yForThisL2 = new std::vector<double>;
      std::vector<double> *zForThisL2 = new std::vector<double>;
      std::vector<int> *stationsForThisL2 = new std::vector<int>;
      int nMuHitsForThisL2 = 0;
      
      // We first start with checking for a valid FTS from the IP: i.e. do we have a valid L2 to seed an L3?  
      // Then we rescale the errorMatrix at the IP.
      
      TrajectoryStateTransform transform; 
      
      // To get our error matrix things, we need a MuonServiceProxy.  Declared here.
      
      MuonServiceProxy *proxy = new MuonServiceProxy(muonServiceParams);
      proxy->update(iSetup);
      
      // And this we need to do because the f.t.s. doesn't take TrackRefs
      const reco::Track *tk = const_cast<const reco::Track *>(&*refL2);
      
      FreeTrajectoryState fts = transform.initialFreeState(*tk,&*proxy->magneticField());
      
      //rescale the error at IP, but only if I have to.  And I have to.
      // matrixForThisL2 is the rescaling matrix to be applied.
      AlgebraicSymMatrix55 matrixForThisL2 = theErrorMatrix->get(GlobalVector(refL2->px(),refL2->py(),refL2->pz())).matrix();
      
      if (theErrorMatrix && theAdjustAtIp){ 
	CurvilinearTrajectoryError oMat = fts.curvilinearError();
	CurvilinearTrajectoryError sfMat = theErrorMatrix->get(fts.momentum());//FIXME with position    
	MuonErrorMatrix::multiply(oMat, sfMat);
	fts = FreeTrajectoryState(fts.parameters(),oMat); 
      }
      
      edm::ESHandle<Propagator> testProp = proxy->propagator(thePropagatorName);
      
      if (testProp.isValid()) {
	
	StateOnTrackerBound onBounds(testProp.product());
	
	TrajectoryStateOnSurface outer = onBounds(fts);
	
	if (theErrorMatrix && !theAdjustAtIp){ 
	  CurvilinearTrajectoryError oMat = outer.curvilinearError();
	  CurvilinearTrajectoryError sfMat = theErrorMatrix->get(outer.globalMomentum());//FIXME with position
	  MuonErrorMatrix::multiply(oMat, sfMat);   
	  outer = TrajectoryStateOnSurface(outer.globalParameters(),
					   oMat,
					   outer.surface(),
					   outer.surfaceSide(),
					   outer.weight());
	}
	
	if (outer.isValid()) {
	  AlgebraicSymMatrix55 matrixForThisL2 = outer.curvilinearError().matrix();
	  for (int i = 0; i < 5; i++) {
	    for (int j = 0; j < 5; j++) {
	      double temp = matrixForThisL2(i,j);
	      (*matrixValuesForThisL2).push_back(temp);
	    }
	  }
	}
	(*muonErrorMatrix).insert(std::make_pair(iL2,*matrixValuesForThisL2));
      }
      matrixValuesForThisL2->clear();
      
      
      
      edm::ESHandle<TransientTrackingRecHitBuilder> muonBuilder;
      std::string muonBuilderName = "MuonRecHitBuilder";
      iSetup.get<TransientRecHitRecord>().get(muonBuilderName,muonBuilder);
      
      for (trackingRecHit_iterator l2Hit = refL2->recHitsBegin(); l2Hit != refL2->recHitsEnd(); ++l2Hit) {
	if ((*l2Hit)->isValid()) {
	  (*idsForThisL2).push_back((*l2Hit)->geographicalId().rawId());
	  (*subidsForThisL2).push_back((*l2Hit)->geographicalId().subdetId());
	  (*detsForThisL2).push_back((*l2Hit)->geographicalId().det());
	  (*statusForThisL2).push_back((*l2Hit)->type());
	  if ((*l2Hit)->geographicalId().det() == 2) { // Muon System
	    nMuHitsForThisL2++;
	    TransientTrackingRecHit::RecHitPointer globL2 = muonBuilder->build(&**l2Hit);
	    (*xForThisL2).push_back(globL2->globalPosition().x());
	    (*yForThisL2).push_back(globL2->globalPosition().y());
	    (*zForThisL2).push_back(globL2->globalPosition().z());
	    if ( (*l2Hit)->geographicalId().subdetId() == 1) { // DT hit
	      const DTChamberId& id = DTChamberId((*l2Hit)->geographicalId());
	      (*stationsForThisL2).push_back(id.station());
	    }
	    if ( (*l2Hit)->geographicalId().subdetId() == 2) { // CSC hit
	      const CSCDetId& id = CSCDetId((*l2Hit)->geographicalId());
	      (*stationsForThisL2).push_back(id.station());
	    }
	    if ( (*l2Hit)->geographicalId().subdetId() == 3) { // RPC hit
	      const RPCDetId& id = RPCDetId((*l2Hit)->geographicalId());
	      (*stationsForThisL2).push_back(id.station());
	    }
	  }
	  else edm::LogError(theCategory)<<"Hits for L2 muon outside muon system.";
	}
      }
      
      (*l2DetIds).insert(std::make_pair(iL2,*idsForThisL2));
      (*l2SubdetIds).insert(std::make_pair(iL2,*subidsForThisL2));
      (*l2Component).insert(std::make_pair(iL2,*detsForThisL2));
      (*l2RecHitsStatus).insert(std::make_pair(iL2,*statusForThisL2));
      (*l2NMuHits).insert(std::make_pair(iL2,nMuHitsForThisL2));
      (*l2MuStationNumber).insert(std::make_pair(iL2,*stationsForThisL2));
      (*l2RecHitsX).insert(std::make_pair(iL2,*xForThisL2));
      (*l2RecHitsY).insert(std::make_pair(iL2,*yForThisL2));
      (*l2RecHitsZ).insert(std::make_pair(iL2,*zForThisL2));
      
      idsForThisL2->clear();
      subidsForThisL2->clear();
      detsForThisL2->clear();
      statusForThisL2->clear();
      stationsForThisL2->clear();
      xForThisL2->clear();
      yForThisL2->clear();
      zForThisL2->clear();
      nMuHitsForThisL2 = 0;
      
      // Determine whether this L2 muon seeds L3
      int seedsL3 = 0;
      for (unsigned int i = 0; i < indexL2SeedingL3->size(); i++) {
	if (iL2 == indexL2SeedingL3->at(i)) seedsL3 = 1;
      }
      (*l2SeedsL3).push_back(seedsL3);
      
      // get the Calorimeter isolation deposits for this L2 muon
      reco::IsoDeposit calDeposit = caloDepositExtractor->deposit(iEvent, iSetup, *refL2);
      // cutting for the L2 muon isolation
      muonisolation::Cuts::CutSpec calo_cuts_here = L2IsoCalCuts(refL2->eta());
      // and deposit for the L2 muon isolation
      double conesize = calo_cuts_here.conesize;
      (*l2CalIsoDeposit).push_back(calDeposit.depositWithin(conesize));
      
      // With the detector-level things filled, time to start doing the associations to sim
      bool associated = false;
      // With the detector-level things filled, time to start doing the associations to sim
      for (reco::RecoToSimCollection::const_iterator findRefL2 = l2RecSimColl.begin(); findRefL2 != l2RecSimColl.end(); ++findRefL2) {
	const edm::RefToBase<reco::Track> & l2RecSimMatch = findRefL2->key;
	int sim_index = 0;
	if (l2RecSimMatch->pt() == refL2->pt()) {
	  associated = true;
	  const std::vector<std::pair<TrackingParticleRef,double> > & tp = findRefL2->val;
	  const TrackingParticleRef & trp = tp.begin()->first;
	  
	  (*l2AssociationVar).push_back(tp.begin()->second);
	  
	  int particle_ID = trp->pdgId();
	  // int myBin = wantMotherBin.GetBinNum(particle_ID);
	  
	  if(abs(particle_ID) == 13){
	    // put in the associated pt,eta,phi
	    (*l2AssociatedSimMuonPt).push_back(trp->pt());
	    (*l2AssociatedSimMuonEta).push_back(trp->eta());
	    (*l2AssociatedSimMuonPhi).push_back(trp->phi());
	    // put in the detIDs for this sim muon
	    std::vector<int> *idsForSimL2 = new std::vector<int>;
	    std::vector<int> *stationsForSim = new std::vector<int>;
	    int simHitCounter = 0;
	    int simMuHitCounter = 0;
	    for (PSimHitContainer::const_iterator l2SimHit = trp->pSimHit_begin(); l2SimHit != trp->pSimHit_end(); ++l2SimHit) {
	      (*idsForSimL2).push_back((*l2SimHit).detUnitId());
	      DetId theDetUnitId(l2SimHit->detUnitId());
	      int detector = theDetUnitId.det();
	      int subdetector = theDetUnitId.subdetId();
	      if (detector == 2) { //Muon system
		simMuHitCounter ++;
		if (subdetector == 1) { //DT
		  const DTChamberId& id = DTChamberId(l2SimHit->detUnitId());
		  (*stationsForSim).push_back(id.station());
		}
		if (subdetector == 2) { //CSC
		  const CSCDetId& id=CSCDetId(l2SimHit->detUnitId());
		  (*stationsForSim).push_back(id.station());
		}
		if (subdetector == 3) { //RPC
		  const RPCDetId& id = RPCDetId(l2SimHit->detUnitId());
		  (*stationsForSim).push_back(id.station());
		}
	      }
	      simHitCounter++;
	    }
	    (*l2AssociatedSimMuonNHits).push_back(simHitCounter);
	    (*l2AssociatedSimMuonDetIds).insert(std::make_pair(sim_index,*idsForSimL2));
	    (*l2AssociatedSimMuonMuStationNumber).insert(std::make_pair(sim_index,*stationsForSim));
	    (*l2AssociatedSimMuonNMuHits).insert(std::make_pair(sim_index,simMuHitCounter));
	    idsForSimL2->clear();
	    stationsForSim->clear();
	    sim_index++;
	    //---------------------- MOTHERHOOD ---------------------------
	    //find the parent of tracking particle
	    for(TrackingParticle::g4t_iterator isimtk = trp->g4Track_begin();isimtk!=trp->g4Track_end();isimtk++)
	      {
		LogDebug(theCategory)<<"I am here 1";
		if(isimtk->type()==13||isimtk->type()==-13)
		  {
		    // This is the sim track for this tracking particle.  Time to put in the parameters
		    FreeTrajectoryState 
		      ftsAtProduction(GlobalPoint(trp->vertex().x(),trp->vertex().y(),trp->vertex().z()),
				      GlobalVector(isimtk->momentum().x(),isimtk->momentum().y(),isimtk->momentum().z()),
				      TrackCharge(trp->charge()),
				      field.product());
		    TSCBLBuilderNoMaterial tscblBuilder;
		    TrajectoryStateClosestToBeamLine tsAtClosestApproach = tscblBuilder(ftsAtProduction,bs);//as in TrackProducerAlgorithm
		    GlobalPoint v1 = tsAtClosestApproach.isValid() ? tsAtClosestApproach.trackStateAtPCA().position() : GlobalPoint();
		    GlobalVector p = tsAtClosestApproach.isValid() ? tsAtClosestApproach.trackStateAtPCA().momentum() : GlobalVector();
		    GlobalPoint v(v1.x()-bs.x0(),v1.y()-bs.y0(),v1.z()-bs.z0());
		    
		    double qoverpSim = tsAtClosestApproach.isValid() ? tsAtClosestApproach.trackStateAtPCA().charge()/p.mag() : 0.;
		    double lambdaSim = M_PI/2-p.theta();
		    double dxySim    = (-v.x()*sin(p.phi())+v.y()*cos(p.phi()));
		    double dzSim     = v.z() - (v.x()*p.x()+v.y()*p.y())/p.perp() * p.z()/p.perp();
		    
		    (*l2AssociatedSimMuonDsz).push_back(dzSim);
		    (*l2AssociatedSimMuonDxy).push_back(dxySim);
		    (*l2AssociatedSimMuonLambda).push_back(lambdaSim);
		    (*l2AssociatedSimMuonQoverP).push_back(qoverpSim);
		    //calculate mother hood
		    MotherSearch mother(&*isimtk, SimTk, SimVtx, hepmc);
		    //FIXME, use reco::Particle mother.mother();
		    //                double pt,eta,phi;
		    //                int parentID;
		    //                int motherBinNumber;
		    if (mother.IsValid()){
		      if (mother.SimIsValid()){
			(*l2ParentID).push_back(mother.Sim_mother->type());
			(*l2MotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Sim_mother->type()));
		      }
		      else {
			(*l2ParentID).push_back(mother.Gen_mother->pdg_id());
			(*l2MotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id()));
		      }
		      //do it once per tracking particle once it succeed
		      break;
		    }
		    else{
		      (*l2ParentID).push_back(0);
		      (*l2MotherBinNumber).push_back(wantMotherBin.GetBinNum(0));
		      edm::LogError(theCategory)<<"tricky muon from TrackingParticle.";
		    }
		  }//sim track is a muon
		else{
		  edm::LogError(theCategory)<<"the sim track attached to the tracking particle is not a muon.";
		  (*l2ParentID).push_back(isimtk->type());
		  (*l2MotherBinNumber).push_back(777);
		}
	      }//loop over SimTrack of tracking particle
	  }//muon associated
	  else{
	    //a reco muon is associated to something else than a muon
	    edm::LogError(theCategory)<<"a reconstructed muon is associated to: "<<particle_ID;
	    (*l2ParentID).push_back(particle_ID);
	    (*l2MotherBinNumber).push_back(-777);
	  }
	}//track has an association
	else{
	  //this track was not associated.
	  edm::LogError(theCategory)<<"a reconstructed muon is not associated.";
	}
      }
      if (associated) (*l2IsAssociated).push_back(1);
      else {
	(*l2IsAssociated).push_back(0);
	(*l2AssociationVar).push_back(-999);
	(*l2AssociatedSimMuonPt).push_back(-999);
	(*l2AssociatedSimMuonEta).push_back(-999);
	(*l2AssociatedSimMuonPhi).push_back(-999);
	(*l2AssociatedSimMuonNHits).push_back(-999);
	(*l2AssociatedSimMuonQoverP).push_back(-999);
	(*l2AssociatedSimMuonLambda).push_back(-999);
	(*l2AssociatedSimMuonDxy).push_back(-999);
	(*l2AssociatedSimMuonDsz).push_back(-999);
	(*l2ParentID).push_back(-666);
	(*l2MotherBinNumber).push_back(-999);
      }
    } else { //loop over l2Muons
      double crap = -999;
      (*l2P).push_back(crap);
      (*l2Px).push_back(crap);
      (*l2Py).push_back(crap);
      (*l2Pz).push_back(crap);
      (*l2Pt).push_back(crap);
      (*l2PtError).push_back(crap);
      (*l2Pt90).push_back(crap);
      (*l2Eta).push_back(crap);
      (*l2EtaError).push_back(crap);
      (*l2Phi).push_back(crap);
      (*l2PhiError).push_back(crap);
      (*l2D0).push_back(crap);
      (*l2D0Error).push_back(crap);
      (*l2NHits).push_back(0);
      (*l2Charge).push_back(crap);
      (*l2Chi2).push_back(crap);
      (*l2Ndof).push_back(crap);
      // Fill in the track fitting parameters (with phi filled in above)
      (*l2Dsz).push_back(crap);
      (*l2DszError).push_back(crap);
      (*l2Dxy).push_back(crap);
      (*l2DxyError).push_back(crap);
      (*l2Lambda).push_back(crap);
      (*l2LambdaError).push_back(crap);
      (*l2Qoverp).push_back(crap);
      (*l2QoverpError).push_back(crap);
      //adam (*l2ErrorMatrix).push_back(NULL);
      
      (*l2DetIds).insert(std::make_pair(iMu,0));    
      (*l2SubdetIds).insert(std::make_pair(iMu,0));
      (*l2Component).insert(std::make_pair(iMu,0));
      (*l2RecHitsStatus).insert(std::make_pair(iMu,0));
      (*l2NMuHits).insert(std::make_pair(iMu,0));
      (*l2MuStationNumber).insert(std::make_pair(iMu,0));
      (*l2RecHitsX).insert(std::make_pair(iMu,0));
      (*l2RecHitsY).insert(std::make_pair(iMu,0));
      (*l2RecHitsZ).insert(std::make_pair(iMu,0));
      
      //(*l2IsoTrackDR).push_back(crap);
      //(*l2IsoTrackDRMinDelPt).push_back(crap);
      //(*l2TrackIsoDeposit).push_back(0);
      
      (*l2IsAssociated).push_back(0);
      (*l2AssociationVar).push_back(crap);
      (*l2AssociationPdgId).push_back(-999);
      (*l2AssociationMyBit).push_back(0);
      (*l2AssociationVtxX).push_back(crap);
      (*l2AssociationVtxY).push_back(crap);
      (*l2AssociationVtxZ).push_back(crap);
      (*l2AssociatedSimMuonPt).push_back(crap);
      (*l2AssociatedSimMuonEta).push_back(crap);
      (*l2AssociatedSimMuonPhi).push_back(crap);
      
      (*l2AssociatedSimMuonNHits).push_back(0);
      (*l2AssociatedSimMuonDetIds).insert(std::make_pair(iMu,0));
      (*l2AssociatedSimMuonMuStationNumber).insert(std::make_pair(iMu,0));
      (*l2AssociatedSimMuonNMuHits).insert(std::make_pair(iMu,0));
      
      //      (*l2IsAssociated).push_back(0);
      //      (*l2AssociationVar).push_back(-999);
      //      (*l2AssociatedSimMuonPt).push_back(-999);
      //      (*l2AssociatedSimMuonEta).push_back(-999);
      //      (*l2AssociatedSimMuonPhi).push_back(-999);
      //      (*l2AssociatedSimMuonNHits).push_back(-999);
      (*l2AssociatedSimMuonQoverP).push_back(-999);
      (*l2AssociatedSimMuonLambda).push_back(-999);
      (*l2AssociatedSimMuonDxy).push_back(-999);
      (*l2AssociatedSimMuonDsz).push_back(-999);
      (*l2ParentID).push_back(-666);
      (*l2MotherBinNumber).push_back(-999);
    }
    
    iMu++;
  } // this brace migrates to swallow GLB STA TK
  
  // Loop over all tracking particles

  int sim_index = 0;
  nSimMuon = 0;
  //  for (TrackingParticleCollection::const_iterator trp = (*TPtracks).begin();
  //       trp != (*TPtracks).end(); ++trp) {

  //  edm::LogInfo("MuonRecoTreeUtility") << "total number of sim particles = " << (*TPtracks).size();
  
  for (unsigned int iSim = 0; iSim != (*TPtracks).size(); iSim++) {
    
    TrackingParticleRef trp(TPtracks, iSim);
    int particle_ID = trp->pdgId();
    if(abs(particle_ID) != 13) continue;
    //    if (abs(particle_ID) != 13) edm::LogInfo("MuonRecoTreeUtility") << "we have a non-muon in the collection.";
    (*simMuonPt).push_back(trp->pt());
    (*simMuonEta).push_back(trp->eta());
    (*simMuonPhi).push_back(trp->phi());
    LogDebug("SpecialBit");
    unsigned int thisBit = getBit(trp);
    (*simMuonMyBit).push_back(thisBit);
    (*simMuonVtxX).push_back(trp->vertex().x());
    (*simMuonVtxY).push_back(trp->vertex().y());
    (*simMuonVtxZ).push_back(trp->vertex().z());
    
    std::vector<int> *idsForSim = new std::vector<int>;
    std::vector<int> *stationsForSim = new std::vector<int>;
    int simHitCounter = 0;
    int simMuHitCounter = 0;
    for (PSimHitContainer::const_iterator simHit = trp->pSimHit_begin(); simHit != trp->pSimHit_end(); ++simHit) {
      (*idsForSim).push_back((*simHit).detUnitId());
      DetId theDetUnitId(simHit->detUnitId());
      int detector = theDetUnitId.det();
      int subdetector = theDetUnitId.subdetId();
      if (detector == 2) { //Muon system
	simMuHitCounter ++;
	if (subdetector == 1) { //DT 
	  const DTChamberId& id = DTChamberId(simHit->detUnitId());
	  (*stationsForSim).push_back(id.station());
	}
	if (subdetector == 2) { //CSC
	  const CSCDetId& id=CSCDetId(simHit->detUnitId());
	  (*stationsForSim).push_back(id.station());
	}
	if (subdetector == 3) { //RPC
	  const RPCDetId& id = RPCDetId(simHit->detUnitId());
	  (*stationsForSim).push_back(id.station());
	}
      }
      simHitCounter++;
    }
    (*simMuonNHits).push_back(simHitCounter);
    (*simMuonDetIds).insert(std::make_pair(sim_index,*idsForSim));
    (*simMuonMuStationNumber).insert(std::make_pair(sim_index,*stationsForSim));
    (*simMuonNMuHits).insert(std::make_pair(sim_index,simMuHitCounter));
    idsForSim->clear();
    stationsForSim->clear();
    
    for(TrackingParticle::g4t_iterator isimtk = trp->g4Track_begin();isimtk!=trp->g4Track_end();isimtk++)  {
      if(isimtk->type()==13||isimtk->type()==-13) {
	// This is the sim track for this tracking particle.  Time to put in the parameters
	FreeTrajectoryState 
	  ftsAtProduction(GlobalPoint(trp->vertex().x(),trp->vertex().y(),trp->vertex().z()),
			  GlobalVector(isimtk->momentum().x(),isimtk->momentum().y(),isimtk->momentum().z()),
			  TrackCharge(trp->charge()),
			  field.product());
	TSCBLBuilderNoMaterial tscblBuilder;
	TrajectoryStateClosestToBeamLine tsAtClosestApproach = tscblBuilder(ftsAtProduction,bs);//as in TrackProducerAlgorithm
	GlobalPoint v1 = tsAtClosestApproach.isValid() ? tsAtClosestApproach.trackStateAtPCA().position() : GlobalPoint();
	GlobalVector p = tsAtClosestApproach.isValid() ? tsAtClosestApproach.trackStateAtPCA().momentum() : GlobalVector();
	GlobalPoint v(v1.x()-bs.x0(),v1.y()-bs.y0(),v1.z()-bs.z0());
	
	double qoverpSim = tsAtClosestApproach.isValid() ? tsAtClosestApproach.trackStateAtPCA().charge()/p.mag() : 0.;
	double lambdaSim = M_PI/2-p.theta();
	double dxySim    = (-v.x()*sin(p.phi())+v.y()*cos(p.phi()));
	double dzSim     = v.z() - (v.x()*p.x()+v.y()*p.y())/p.perp() * p.z()/p.perp();
	
	(*simMuonDsz).push_back(dzSim);
	(*simMuonDxy).push_back(dxySim);
	(*simMuonLambda).push_back(lambdaSim);
	(*simMuonQoverP).push_back(qoverpSim);
	//calculate mother hood
	MotherSearch mother(&*isimtk, SimTk, SimVtx, hepmc);
	if (mother.IsValid()){
	  if (mother.SimIsValid()){
	    (*simMuonParentID).push_back(mother.Sim_mother->type());
	    (*simMuonMotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Sim_mother->type()));
	  } // simIsValid
	  else {
	    (*simMuonParentID).push_back(mother.Gen_mother->pdg_id());
	    (*simMuonMotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id()));
	  } // gen used otherwise
	} // motherIsValid
	else {
	  (*simMuonParentID).push_back(0);
	  (*simMuonMotherBinNumber).push_back(wantMotherBin.GetBinNum(0));
	}
	//do it once per tracking particle once it succeed
	break;
      } // sim track is muon
    } //loop over g4_iterator
      // SimToReco associations. First the TrackingParticleRef
      // First look to see if there's a SimToRec match at L3
      /* Temporary commentary here:
	 
      We do an on-the-fly SimToReco association.  There's a small number of events where the result 
      is double-counting: one sim can be associated to 2 reco, if there's a low-pt reco within Delta R 
      of the sim.  The size of this effect is about 1/500 for L2, less for L3.
      
      The MultiTrackValidator uses an AssociatorMap to get around this problem.  I suppose the thing to do
      is to look at what it does, and what I can do to get the same results.  To the LXR, then.
      
      */
    
    std::vector<std::pair<RefToBase<reco::Track>, double> > rt;

    if(l3SimRecColl.find(trp) != l3SimRecColl.end()){
      rt = (std::vector<std::pair<RefToBase<reco::Track>, double> >) l3SimRecColl[trp];
      if (rt.size()!=0) { // Association to L3 successful
	(*simToL3Associated).push_back(1);
	(*simToL3AssociationVar).push_back(rt.begin()->second);
	int iL3 = 0;
	for(View<Muon>::const_iterator iMuon = muonColl.begin();
	    iMuon != muonColl.end(); ++iMuon) {
	  if (iMuon->combinedMuon().isAvailable() && iMuon->combinedMuon()->pt() == (rt.begin()->first->pt() ) ) {
	    (*simToL3RecoIndex).push_back(iL3);
	  }
	  iL3++;
	}
      }
      else { // Something went wrong
	edm::LogInfo("MuonRecoTreeUtility")<<"rt.size = 0, but l3SimRec finds trp";
	(*simToL3Associated).push_back(0);
	(*simToL3AssociationVar).push_back(-999);
	(*simToL3RecoIndex).push_back(-999);
      }
    }
    else { // Association to L3 unsuccessful
      (*simToL3Associated).push_back(0);
      (*simToL3AssociationVar).push_back(-999);
      (*simToL3RecoIndex).push_back(-999);
    }

    if(tkSimRecColl.find(trp) != tkSimRecColl.end()){
      rt = (std::vector<std::pair<RefToBase<reco::Track>, double> >) tkSimRecColl[trp];
      if (rt.size()!=0) { // Association to TK successful
        (*simToTkAssociated).push_back(1);
        (*simToTkAssociationVar).push_back(rt.begin()->second);
	int iTk = 0;
	for(View<Muon>::const_iterator iMuon = muonColl.begin();
	    iMuon != muonColl.end(); ++iMuon) {
	  if (iMuon->innerTrack().isAvailable() && iMuon->innerTrack()->pt() == (rt.begin()->first->pt() ) ) {
            (*simToTkRecoIndex).push_back(iTk);
          }
	  iTk++;
        }
      }
      else { // Something went wrong
	edm::LogInfo("MuonRecoTreeUtility")<<"rt.size = 0, but l3SimRec finds trp";
        (*simToTkAssociated).push_back(0);
        (*simToTkAssociationVar).push_back(-999);
        (*simToTkRecoIndex).push_back(-999);
      }
    }
    else { // Association to L3 unsuccessful
      (*simToTkAssociated).push_back(0);
      (*simToTkAssociationVar).push_back(-999);
      (*simToTkRecoIndex).push_back(-999);
    }

    if(l2SimRecColl.find(trp) != l2SimRecColl.end()){
      rt = (std::vector<std::pair<RefToBase<reco::Track>, double> >) l2SimRecColl[trp];
      if (rt.size()!=0) { // Association to L2 successful
	(*simToL2Associated).push_back(1);
	(*simToL2AssociationVar).push_back(rt.begin()->second);
	int iL2 = 0;
	for(View<Muon>::const_iterator iMuon = muonColl.begin();
	    iMuon != muonColl.end(); ++iMuon) {
	  if (iMuon->outerTrack().isAvailable() && iMuon->outerTrack()->pt() == (rt.begin()->first->pt() ) ) {
	    (*simToL2RecoIndex).push_back(iL2);
	  }
	}
      }
      else { // Something went wrong
	edm::LogInfo("MuonRecoTreeUtility")<<"rt.size = 0, but l2SimRec finds trp";
	(*simToL2Associated).push_back(0);
	(*simToL2AssociationVar).push_back(-999);
	(*simToL2RecoIndex).push_back(-999);
      }
    }
    else { // Association to L2 unsuccessful
      (*simToL2Associated).push_back(0);
      (*simToL2AssociationVar).push_back(-999);
      (*simToL2RecoIndex).push_back(-999);
    }
    
    sim_index++;
    
    //    } //trackingParticle is a muon
  }//loop over all trackingparticles
  
  nSimMuon = sim_index;

  MuTrigData->Fill();
  MuTrigMC->Fill();

  /*
  triggerDecisions->clear(); 
  triggerNames->clear();
  */
  /*
  muonDigiModuleTimes->clear();
  muonLocalRecModuleTimes->clear();
  muonL2RecModuleTimes->clear();
  muonL3RecModuleTimes->clear();
  muonL2IsoModuleTimes->clear();
  muonL3IsoModuleTimes->clear();
  trackerDigiModuleTimes->clear();
  trackerRecModuleTimes->clear();
  caloDigiModuleTimes->clear();
  caloRecModuleTimes->clear();
  */

  muAllGlobalMuons->clear();
  muAllStandAloneMuons ->clear();
  muAllTrackerMuons ->clear();
  muTrackerMuonArbitrated ->clear();
  muAllArbitrated ->clear();
  muGlobalMuonPromptTight ->clear();
  muTMLastStationLoose ->clear();
  muTMLastStationTight ->clear();
  muTM2DCompatibilityLoose ->clear();
  muTM2DCompatibilityTight ->clear();
  muTMOneStationLoose ->clear();
  muTMOneStationTight ->clear();
  muTMLastStationOptimizedLowPtLoose ->clear();
  muTMLastStationOptimizedLowPtTight ->clear();
  muGMTkChiCompatibility ->clear();
  muGMStaChiCompatibility ->clear();
  muGMTkKinkTight ->clear();
  muCaloCompatibility ->clear();
  muSegmentCompatibility ->clear();
  muTrkKink ->clear();
  muGlbKink ->clear();
  muTrkRelChi2 ->clear();
  muStaRelChi2 ->clear();
  //  muIso03Valid->clear();
  //  muIso03sumPt->clear();
  //  muIso03emEt->clear();
  //  muIso03hadEt->clear();
  //  muIso03nTracks->clear();
  //  muIso03trackerVetoPt->clear();
  muNumberOfChambers->clear();
  muNumberOfMatches->clear();
  muStationMask->clear();

  muNCSCSeg->clear();
  muNDTSeg->clear();
  muNRPCSeg->clear();

  l3P->clear();
  l3Px->clear();
  l3Py->clear();
  l3Pz->clear();
  l3Pt->clear();
  l3PtError->clear();
  l3Pt90->clear();
  l3Eta->clear();
  l3EtaError->clear();
  l3Phi->clear();
  l3PhiError->clear();
  l3D0->clear();
  l3D0Error->clear();
  l3NHits->clear();
  l3Charge->clear();
  l3Chi2->clear();
  l3Ndof->clear();
  l3DetIds->clear();
  l3SubdetIds->clear();
  l3Component->clear();
  l3NMuHits->clear();
  l3MuStationNumber->clear();
  l3RecHitsStatus->clear();
  l3RecHitsX->clear();
  l3RecHitsY->clear();
  l3RecHitsZ->clear();
  l3RecHitsXTM->clear();
  l3RecHitsYTM->clear();
  l3RecHitsZTM->clear();
  l3RecHitsXTSOS->clear();
  l3RecHitsYTSOS->clear();
  l3RecHitsZTSOS->clear();
  l3RecHitsPhiTM->clear();
  l3RecHitsErrorTM->clear();
  l3RecHitsPhiTSOS->clear();

  l3CalIsoDeposit->clear();
  l3TrackIsoDeposit->clear();
  l3IsoTrackDR->clear();
  l3IsoTrackDRMinDelPt->clear();

  l3Dsz->clear();
  l3DszError->clear();
  l3Dxy->clear();
  l3DxyError->clear();
  l3Lambda->clear();
  l3LambdaError->clear();
  l3Qoverp->clear();
  l3QoverpError->clear();
  l3ErrorMatrix->clear();

  l2SeedsL3->clear();
  indexL2SeedingL3->clear();
  indexL3SeededFromL2->clear();

  muonErrorMatrix->clear();

  l3IsAssociated->clear();
  l3ParentID->clear();
  l3MotherBinNumber->clear();
  l3AssociationVar->clear();
  l3AssociationPdgId->clear();
  l3AssociationMyBit->clear();
  l3AssociationVtxX->clear();
  l3AssociationVtxY->clear();
  l3AssociationVtxZ->clear();
  l3AssociatedSimMuonIndex->clear();
  l3AssociatedSimMuonPt->clear();
  l3AssociatedSimMuonEta->clear();
  l3AssociatedSimMuonPhi->clear();
  l3AssociatedSimMuonNHits->clear();
  l3AssociatedSimMuonDetIds->clear();
  l3AssociatedSimMuonNMuHits->clear();
  l3AssociatedSimMuonMuStationNumber->clear();

  l3AssociatedSimMuonDsz->clear();
  l3AssociatedSimMuonDxy->clear();
  l3AssociatedSimMuonLambda->clear();
  l3AssociatedSimMuonQoverP->clear();

  tkTrackP->clear();
  tkTrackPx->clear();
  tkTrackPy->clear();
  tkTrackPz->clear();
  tkTrackPt->clear();
  tkTrackPtError->clear();
  tkTrackEta->clear();
  tkTrackEtaError->clear();
  tkTrackPhi->clear();
  tkTrackPhiError->clear();
  tkTrackD0->clear();
  tkTrackD0Error->clear();
  tkTrackNHits->clear();
  tkTrackCharge->clear();
  tkTrackChi2->clear();
  tkTrackNdof->clear();
  tkTrackDetIds->clear();
  tkTrackSubdetIds->clear();
  tkTrackRecHitsStatus->clear();
  tkTrackRecHitsX->clear();
  tkTrackRecHitsY->clear();
  tkTrackRecHitsZ->clear();

  tkTrackDsz->clear();
  tkTrackDszError->clear();
  tkTrackDxy->clear();
  tkTrackDxyError->clear();
  tkTrackLambda->clear();
  tkTrackLambdaError->clear();
  tkTrackQoverp->clear();
  tkTrackQoverpError->clear();
  tkTrackErrorMatrix->clear();

  tkTrackIsAssociated->clear();
  tkTrackParentID->clear();
  tkTrackMotherBinNumber->clear();
  tkTrackAssociationVar->clear();
  tkTrackAssociationPdgId->clear();
  tkTrackAssociationMyBit->clear();
  tkTrackAssociationVtxX->clear();
  tkTrackAssociationVtxY->clear();
  tkTrackAssociationVtxZ->clear();
  tkTrackAssociatedSimMuonPt->clear();
  tkTrackAssociatedSimMuonEta->clear();
  tkTrackAssociatedSimMuonPhi->clear();
  tkTrackAssociatedSimMuonNHits->clear();
  tkTrackAssociatedSimMuonDetIds->clear();

  tkTrackAssociatedSimMuonDsz->clear();
  tkTrackAssociatedSimMuonDxy->clear();
  tkTrackAssociatedSimMuonLambda->clear();
  tkTrackAssociatedSimMuonQoverP->clear();

  l2P->clear();
  l2Px->clear();
  l2Py->clear();
  l2Pz->clear();
  l2Pt->clear();
  l2PtError->clear();
  l2Pt90->clear();
  l2Eta->clear();
  l2EtaError->clear();
  l2Phi->clear();
  l2PhiError->clear();
  l2D0->clear();
  l2D0Error->clear();
  l2NHits->clear();
  l2Charge->clear();
  l2Chi2->clear();
  l2Ndof->clear();
  l2DetIds->clear();
  l2SubdetIds->clear();
  l2Component->clear();
  l2NMuHits->clear();
  l2MuStationNumber->clear();
  l2RecHitsStatus->clear();
  l2RecHitsX->clear();
  l2RecHitsY->clear();
  l2RecHitsZ->clear();

  l2CalIsoDeposit->clear();

  l2Dsz->clear();
  l2DszError->clear();
  l2Dxy->clear();
  l2DxyError->clear();
  l2Lambda->clear();
  l2LambdaError->clear();
  l2Qoverp->clear();
  l2QoverpError->clear();
  l2ErrorMatrix->clear();

  l2IsAssociated->clear();
  l2ParentID->clear();
  l2MotherBinNumber->clear();
  l2AssociationVar->clear();
  l2AssociationPdgId->clear();
  l2AssociationMyBit->clear();
  l2AssociationVtxX->clear();
  l2AssociationVtxY->clear();
  l2AssociationVtxZ->clear();
  l2AssociatedSimMuonPt->clear();
  l2AssociatedSimMuonEta->clear();
  l2AssociatedSimMuonPhi->clear();
  l2AssociatedSimMuonNHits->clear();
  l2AssociatedSimMuonDetIds->clear();
  l2AssociatedSimMuonNMuHits->clear();
  l2AssociatedSimMuonMuStationNumber->clear();

  l2AssociatedSimMuonDsz->clear();
  l2AssociatedSimMuonDxy->clear();
  l2AssociatedSimMuonLambda->clear();
  l2AssociatedSimMuonQoverP->clear();

  simMuonParentID->clear();
  simMuonMotherBinNumber->clear();
  simMuonPt->clear();
  simMuonEta->clear();
  simMuonPhi->clear();
  simMuonMyBit->clear();
  simMuonVtxX->clear();
  simMuonVtxY->clear();
  simMuonVtxZ->clear();
  simMuonNHits->clear();
  simMuonDetIds->clear();
  simMuonNMuHits->clear();
  simMuonMuStationNumber->clear();

  simMuonDsz->clear();
  simMuonDxy->clear();
  simMuonLambda->clear();
  simMuonQoverP->clear();

  simToL3Associated->clear();
  simToL3AssociationVar->clear();
  simToL3RecoIndex->clear();
  simToTkAssociated->clear();
  simToTkAssociationVar->clear();
  simToTkRecoIndex->clear();
  simToL2Associated->clear();
  simToL2AssociationVar->clear();
  simToL2RecoIndex->clear();


}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonRecoTreeUtility::beginJob(const edm::EventSetup&)
{

  theFile = new TFile(outputFileName.c_str(),"recreate");
  theFile->cd();
  
  MuTrigData = new TTree("MuTrigData","MuTrigData");
  MuTrigMC = new TTree("MuTrigMC","MuTrigMC");

  // MuTrigData branches

  //Event-level information
  MuTrigData->Branch("RunNumber",&RunNumber,"RunNumber/I");
  MuTrigData->Branch("EventNumber",&EventNumber,"EventNumber/I");
  /*
  // Execution times
  MuTrigData->Branch("totalMuonHLTTime",&totalMuonHLTTime,"totalMuonHLTTime/D");
  MuTrigData->Branch("muonDigiModuleTimes",&muonDigiModuleTimes);  
  MuTrigData->Branch("muonLocalRecModuleTimes",&muonLocalRecModuleTimes);
  MuTrigData->Branch("muonL2RecModuleTimes",&muonL2RecModuleTimes);
  MuTrigData->Branch("muonL3RecModuleTimes",&muonL3RecModuleTimes);
  MuTrigData->Branch("muonL2IsoModuleTimes",&muonL2IsoModuleTimes);
  MuTrigData->Branch("muonL3IsoModuleTimes",&muonL3IsoModuleTimes);
  MuTrigData->Branch("trackerDigiModuleTimes",&trackerDigiModuleTimes);
  MuTrigData->Branch("trackerRecModuleTimes",&trackerRecModuleTimes);
  MuTrigData->Branch("caloDigiModuleTimes",&caloDigiModuleTimes);
  MuTrigData->Branch("caloRecModuleTimes",&caloRecModuleTimes);
  */
  /*
  //Trigger information
  MuTrigData->Branch("l1SingleMuNonIsoTriggered",&l1SingleMuNonIsoTriggered,"l1SingleMuNonIsoTriggered/I");
  MuTrigData->Branch("l2SingleMuNonIsoTriggered",&l2SingleMuNonIsoTriggered,"l2SingleMuNonIsoTriggered/I");
  MuTrigData->Branch("l3SingleMuNonIsoTriggered",&l3SingleMuNonIsoTriggered,"l3SingleMuNonIsoTriggered/I");
  MuTrigData->Branch("l1SingleMuIsoTriggered",&l1SingleMuIsoTriggered,"l1SingleMuIsoTriggered/I");
  MuTrigData->Branch("l2SingleMuIsoPreTriggered",&l2SingleMuIsoPreTriggered,"l2SingleMuIsoPreTriggered/I");
  MuTrigData->Branch("l2SingleMuIsoTriggered",&l2SingleMuIsoTriggered,"l2SingleMuIsoTriggered/I");
  MuTrigData->Branch("l3SingleMuIsoPreTriggered",&l3SingleMuIsoPreTriggered,"l3SingleMuIsoPreTriggered/I");
  MuTrigData->Branch("l3SingleMuIsoTriggered",&l3SingleMuIsoTriggered,"l3SingleMuIsoTriggered/I");
  MuTrigData->Branch("l1DiMuNonIsoTriggered",&l1DiMuNonIsoTriggered,"l1DiMuNonIsoTriggered/I");
  MuTrigData->Branch("l2DiMuNonIsoTriggered",&l2DiMuNonIsoTriggered,"l2DiMuNonIsoTriggered/I");
  MuTrigData->Branch("l3DiMuNonIsoTriggered",&l3DiMuNonIsoTriggered,"l3DiMuNonIsoTriggered/I");
  MuTrigData->Branch("l1DiMuIsoTriggered",&l1DiMuIsoTriggered,"l1DiMuIsoTriggered/I");
  MuTrigData->Branch("l2DiMuIsoPreTriggered",&l2DiMuIsoPreTriggered,"l2DiMuIsoPreTriggered/I");
  MuTrigData->Branch("l2DiMuIsoTriggered",&l2DiMuIsoTriggered,"l2DiMuIsoTriggered/I");
  MuTrigData->Branch("l3DiMuIsoPreTriggered",&l3DiMuIsoPreTriggered,"l3DiMuIsoPreTriggered/I");
  MuTrigData->Branch("l3DiMuIsoTriggered",&l3DiMuIsoTriggered,"l3DiMuIsoTriggered/I"); 
  MuTrigData->Branch("triggerDecisions",&triggerDecisions);
  MuTrigData->Branch("triggerNames",&triggerNames);
  */
  // number of muons at each level
  MuTrigData->Branch("nMu",&nMu,"nMu/I");
  MuTrigData->Branch("nL2",&nL2,"nL2/I");
  MuTrigData->Branch("nL3",&nL3,"nL3/I");
  MuTrigData->Branch("nTkTracks",&nTkTracks,"nTkTracks/I");
  MuTrigData->Branch("nL3Cands",&nL3Cands,"nL3Cands/I");
  MuTrigData->Branch("muAllGlobalMuons",&muAllGlobalMuons);
  MuTrigData->Branch("muAllStandAloneMuons",&muAllStandAloneMuons);
  MuTrigData->Branch("muAllTrackerMuons",&muAllTrackerMuons);
  MuTrigData->Branch("muTrackerMuonArbitrated",&muTrackerMuonArbitrated);
  MuTrigData->Branch("muAllArbitrated",&muAllArbitrated);
  MuTrigData->Branch("muGlobalMuonPromptTight",&muGlobalMuonPromptTight);
  MuTrigData->Branch("muTMLastStationLoose",&muTMLastStationLoose);
  MuTrigData->Branch("muTMLastStationTight",&muTMLastStationTight);
  MuTrigData->Branch("muTM2DCompatibilityLoose",&muTM2DCompatibilityLoose);
  MuTrigData->Branch("muTM2DCompatibilityTight",&muTM2DCompatibilityTight);
  MuTrigData->Branch("muTMOneStationLoose",&muTMOneStationLoose);
  MuTrigData->Branch("muTMOneStationTight",&muTMOneStationTight);
  MuTrigData->Branch("muTMLastStationOptimizedLowPtLoose",&muTMLastStationOptimizedLowPtLoose);
  MuTrigData->Branch("muTMLastStationOptimizedLowPtTight",&muTMLastStationOptimizedLowPtTight);
  MuTrigData->Branch("muGMTkChiCompatibility",&muGMTkChiCompatibility);
  MuTrigData->Branch("muGMStaChiCompatibility",&muGMStaChiCompatibility);
  MuTrigData->Branch("muGMTkKinkTight",&muGMTkKinkTight);
  MuTrigData->Branch("muCaloCompatibility",&muCaloCompatibility);
  MuTrigData->Branch("muSegmentCompatibility",&muSegmentCompatibility);
  MuTrigData->Branch("muTrkKink",&muTrkKink);
  MuTrigData->Branch("muGlbKink",&muGlbKink);
  MuTrigData->Branch("muTrkRelChi2",&muTrkRelChi2);
  MuTrigData->Branch("muStaRelChi2",&muStaRelChi2);
  //  MuTrigData->Branch("muIso03Valid",&muIso03Valid);
  //  MuTrigData->Branch("muIso03sumPt",&muIso03sumPt);
  //  MuTrigData->Branch("muIso03emEt",&muIso03emEt);
  //  MuTrigData->Branch("muIso03hadEt",&muIso03hadEt);
  //  MuTrigData->Branch("muIso03nTracks",&muIso03nTracks);
  //  MuTrigData->Branch("muIso03trackerVetoPt",&muIso03trackerVetoPt);
  MuTrigData->Branch("muNumberOfChambers",&muNumberOfChambers);
  MuTrigData->Branch("muNumberOfMatches",&muNumberOfMatches);
  MuTrigData->Branch("muStationMask",&muStationMask);
  MuTrigData->Branch("muNCSCSeg",&muNCSCSeg);
  MuTrigData->Branch("muNDTSeg",&muNDTSeg);
  MuTrigData->Branch("muNRPCSeg",&muNRPCSeg);
  //MuTrigData->Branch("nL3Seeds",&nL3Seeds,"nL3Seeds/I");
  // L3 muon information: the basics
  MuTrigData->Branch("l3P",&l3P);
  MuTrigData->Branch("l3Px",&l3Px);
  MuTrigData->Branch("l3Py",&l3Py);
  MuTrigData->Branch("l3Pz",&l3Pz);
  MuTrigData->Branch("l3Pt",&l3Pt);
  MuTrigData->Branch("l3PtError",&l3PtError);
  MuTrigData->Branch("l3Pt90",&l3Pt90);
  MuTrigData->Branch("l3Eta",&l3Eta);
  MuTrigData->Branch("l3EtaError",&l3EtaError);
  MuTrigData->Branch("l3Phi",&l3Phi);
  MuTrigData->Branch("l3PhiError",&l3PhiError);
  MuTrigData->Branch("l3D0",&l3D0);
  MuTrigData->Branch("l3D0Error",&l3D0Error);
  MuTrigData->Branch("l3NHits",&l3NHits);
  MuTrigData->Branch("l3Charge",&l3Charge);
  MuTrigData->Branch("l3Chi2",&l3Chi2);
  MuTrigData->Branch("l3Ndof",&l3Ndof);
  // L3 Muon Isolation quantities
  MuTrigData->Branch("l3CalIsoDeposit",&l3CalIsoDeposit);
  MuTrigData->Branch("l3TrackIsoDeposit",&l3TrackIsoDeposit);
  MuTrigData->Branch("l3IsoTrackDR",&l3IsoTrackDR);
  MuTrigData->Branch("l3IsoTrackDRMinDelPt",&l3IsoTrackDRMinDelPt);
  // L3 Muon Track fitting parameters (with phi already declared above)
  MuTrigData->Branch("l3Dsz",&l3Dsz);
  MuTrigData->Branch("l3DszError",&l3DszError);
  MuTrigData->Branch("l3Dxy",&l3Dxy);
  MuTrigData->Branch("l3DxyError",&l3DxyError);
  MuTrigData->Branch("l3Lambda",&l3Lambda);
  MuTrigData->Branch("l3LambdaError",&l3LambdaError);
  MuTrigData->Branch("l3Qoverp",&l3Qoverp);
  MuTrigData->Branch("l3QoverpError",&l3QoverpError);
  MuTrigData->Branch("l3ErrorMatrix",&l3ErrorMatrix);
  MuTrigData->Branch("l3DetIds",&l3DetIds);
  MuTrigData->Branch("l3SubdetIds",&l3SubdetIds);
  MuTrigData->Branch("l3Component",&l3Component);
  MuTrigData->Branch("l3NMuHits",&l3NMuHits);
  MuTrigData->Branch("l3MuStationNumber",&l3MuStationNumber);
  MuTrigData->Branch("l3RecHitsStatus",&l3RecHitsStatus);
  MuTrigData->Branch("l3RecHitsX",&l3RecHitsX);
  MuTrigData->Branch("l3RecHitsY",&l3RecHitsY);
  MuTrigData->Branch("l3RecHitsZ",&l3RecHitsZ);
  MuTrigData->Branch("l3RecHitsXTM",&l3RecHitsXTM);
  MuTrigData->Branch("l3RecHitsYTM",&l3RecHitsYTM);
  MuTrigData->Branch("l3RecHitsZTM",&l3RecHitsZTM);
  MuTrigData->Branch("l3RecHitsXTSOS",&l3RecHitsXTSOS);
  MuTrigData->Branch("l3RecHitsYTSOS",&l3RecHitsYTSOS);
  MuTrigData->Branch("l3RecHitsZTSOS",&l3RecHitsZTSOS);
  MuTrigData->Branch("l3RecHitsPhiTM",&l3RecHitsPhiTM);
  MuTrigData->Branch("l3RecHitsErrorTM",&l3RecHitsErrorTM);
  MuTrigData->Branch("l3RecHitsPhiTSOS",&l3RecHitsPhiTSOS);
  // Indices for L3<->L2
  MuTrigData->Branch("indexL2SeedingL3",&indexL2SeedingL3);
  MuTrigData->Branch("indexL3SeededFromL2",&indexL3SeededFromL2);
  MuTrigData->Branch("l2SeedsL3",&l2SeedsL3);

  // The muon error matrix
  MuTrigData->Branch("muonErrorMatrix",&muonErrorMatrix);

  // information from hltL3TkTracksFromL2
  MuTrigData->Branch("tkTrackP",&tkTrackP);
  MuTrigData->Branch("tkTrackPx",&tkTrackPx);
  MuTrigData->Branch("tkTrackPy",&tkTrackPy);
  MuTrigData->Branch("tkTrackPz",&tkTrackPz);
  MuTrigData->Branch("tkTrackPt",&tkTrackPt);
  MuTrigData->Branch("tkTrackPtError",&tkTrackPtError);
  MuTrigData->Branch("tkTrackEta",&tkTrackEta);
  MuTrigData->Branch("tkTrackEtaError",&tkTrackEtaError);
  MuTrigData->Branch("tkTrackPhi",&tkTrackPhi);
  MuTrigData->Branch("tkTrackPhiError",&tkTrackPhiError);
  MuTrigData->Branch("tkTrackD0",&tkTrackD0);
  MuTrigData->Branch("tkTrackD0Error",&tkTrackD0Error);
  MuTrigData->Branch("tkTrackNHits",&tkTrackNHits);
  MuTrigData->Branch("tkTrackCharge",&tkTrackCharge);
  MuTrigData->Branch("tkTrackChi2",&tkTrackChi2);
  MuTrigData->Branch("tkTrackNdof",&tkTrackNdof);
  // L3TRACK Muon Track fitting parameters (with phi already declared above)
  MuTrigData->Branch("tkTrackDsz",&tkTrackDsz);
  MuTrigData->Branch("tkTrackDszError",&tkTrackDszError);
  MuTrigData->Branch("tkTrackDxy",&tkTrackDxy);
  MuTrigData->Branch("tkTrackDxyError",&tkTrackDxyError);
  MuTrigData->Branch("tkTrackLambda",&tkTrackLambda);
  MuTrigData->Branch("tkTrackLambdaError",&tkTrackLambdaError);
  MuTrigData->Branch("tkTrackQoverp",&tkTrackQoverp);
  MuTrigData->Branch("tkTrackQoverpError",&tkTrackQoverpError);
  MuTrigData->Branch("tkTrackErrorMatrix",&tkTrackErrorMatrix);
  MuTrigData->Branch("tkTrackDetIds",&tkTrackDetIds);
  MuTrigData->Branch("tkTrackSubdetIds",&tkTrackSubdetIds);
  MuTrigData->Branch("tkTrackRecHitsStatus",&tkTrackRecHitsStatus);
  MuTrigData->Branch("tkTrackRecHitsX",&tkTrackRecHitsX);
  MuTrigData->Branch("tkTrackRecHitsY",&tkTrackRecHitsY);
  MuTrigData->Branch("tkTrackRecHitsZ",&tkTrackRecHitsZ);

  // L2 muon information: the basics
  MuTrigData->Branch("l2P",&l2P);
  MuTrigData->Branch("l2Px",&l2Px);
  MuTrigData->Branch("l2Py",&l2Py);
  MuTrigData->Branch("l2Pz",&l2Pz);
  MuTrigData->Branch("l2Pt",&l2Pt);
  MuTrigData->Branch("l2PtError",&l2PtError);
  MuTrigData->Branch("l2Pt90",&l2Pt90);
  MuTrigData->Branch("l2Eta",&l2Eta);
  MuTrigData->Branch("l2EtaError",&l2EtaError);
  MuTrigData->Branch("l2Phi",&l2Phi);
  MuTrigData->Branch("l2PhiError",&l2PhiError);
  MuTrigData->Branch("l2D0",&l2D0);
  MuTrigData->Branch("l2D0Error",&l2D0Error);
  MuTrigData->Branch("l2NHits",&l2NHits);
  MuTrigData->Branch("l2Charge",&l2Charge);
  MuTrigData->Branch("l2Chi2",&l2Chi2);
  MuTrigData->Branch("l2Ndof",&l2Ndof);
  // L2 Muon Isolation quantities
  MuTrigData->Branch("l2CalIsoDeposit",&l2CalIsoDeposit);
  // L2 Muon Track fitting parameters (with phi already declared above)
  MuTrigData->Branch("l2Dsz",&l2Dsz);
  MuTrigData->Branch("l2DszError",&l2DszError);
  MuTrigData->Branch("l2Dxy",&l2Dxy);
  MuTrigData->Branch("l2DxyError",&l2DxyError);
  MuTrigData->Branch("l2Lambda",&l2Lambda);
  MuTrigData->Branch("l2LambdaError",&l2LambdaError);
  MuTrigData->Branch("l2Qoverp",&l2Qoverp);
  MuTrigData->Branch("l2QoverpError",&l2QoverpError);
  MuTrigData->Branch("l2ErrorMatrix",&l2ErrorMatrix);
  MuTrigData->Branch("l2DetIds",&l2DetIds);
  MuTrigData->Branch("l2SubdetIds",&l2SubdetIds);
  MuTrigData->Branch("l2Component",&l2Component);
  MuTrigData->Branch("l2NMuHits",&l2NMuHits);
  MuTrigData->Branch("l2MuStationNumber",&l2MuStationNumber);
  MuTrigData->Branch("l2RecHitsStatus",&l2RecHitsStatus);
  MuTrigData->Branch("l2RecHitsX",&l2RecHitsX);
  MuTrigData->Branch("l2RecHitsY",&l2RecHitsY);
  MuTrigData->Branch("l2RecHitsZ",&l2RecHitsZ);

  // MuTrigMC branches

  //Event-level information
  MuTrigMC->Branch("RunNumber",&RunNumber,"RunNumber/I");
  MuTrigMC->Branch("EventNumber",&EventNumber,"EventNumber/I");
  /*
  // Execution times
  MuTrigMC->Branch("totalMuonHLTTime",&totalMuonHLTTime,"totalMuonHLTTime/D");
  MuTrigMC->Branch("muonDigiModuleTimes",&muonDigiModuleTimes);
  MuTrigMC->Branch("muonLocalRecModuleTimes",&muonLocalRecModuleTimes);
  MuTrigMC->Branch("muonL2RecModuleTimes",&muonL2RecModuleTimes);
  MuTrigMC->Branch("muonL3RecModuleTimes",&muonL3RecModuleTimes);
  MuTrigMC->Branch("muonL2IsoModuleTimes",&muonL2IsoModuleTimes);
  MuTrigMC->Branch("muonL3IsoModuleTimes",&muonL3IsoModuleTimes);
  MuTrigMC->Branch("trackerDigiModuleTimes",&trackerDigiModuleTimes);
  MuTrigMC->Branch("trackerRecModuleTimes",&trackerRecModuleTimes);
  MuTrigMC->Branch("caloDigiModuleTimes",&caloDigiModuleTimes);
  MuTrigMC->Branch("caloRecModuleTimes",&caloRecModuleTimes);
  */
  /*
  //Trigger information
  MuTrigMC->Branch("l1SingleMuNonIsoTriggered",&l1SingleMuNonIsoTriggered,"l1SingleMuNonIsoTriggered/I");
  MuTrigMC->Branch("l2SingleMuNonIsoTriggered",&l2SingleMuNonIsoTriggered,"l2SingleMuNonIsoTriggered/I");
  MuTrigMC->Branch("l3SingleMuNonIsoTriggered",&l3SingleMuNonIsoTriggered,"l3SingleMuNonIsoTriggered/I");
  MuTrigMC->Branch("l1SingleMuIsoTriggered",&l1SingleMuIsoTriggered,"l1SingleMuIsoTriggered/I");
  MuTrigMC->Branch("l2SingleMuIsoPreTriggered",&l2SingleMuIsoPreTriggered,"l2SingleMuIsoPreTriggered/I");
  MuTrigMC->Branch("l2SingleMuIsoTriggered",&l2SingleMuIsoTriggered,"l2SingleMuIsoTriggered/I");
  MuTrigMC->Branch("l3SingleMuIsoPreTriggered",&l3SingleMuIsoPreTriggered,"l3SingleMuIsoPreTriggered/I");
  MuTrigMC->Branch("l3SingleMuIsoTriggered",&l3SingleMuIsoTriggered,"l3SingleMuIsoTriggered/I");
  MuTrigMC->Branch("l1DiMuNonIsoTriggered",&l1DiMuNonIsoTriggered,"l1DiMuNonIsoTriggered/I");
  MuTrigMC->Branch("l2DiMuNonIsoTriggered",&l2DiMuNonIsoTriggered,"l2DiMuNonIsoTriggered/I");
  MuTrigMC->Branch("l3DiMuNonIsoTriggered",&l3DiMuNonIsoTriggered,"l3DiMuNonIsoTriggered/I");
  MuTrigMC->Branch("l1DiMuIsoTriggered",&l1DiMuIsoTriggered,"l1DiMuIsoTriggered/I");
  MuTrigMC->Branch("l2DiMuIsoPreTriggered",&l2DiMuIsoPreTriggered,"l2DiMuIsoPreTriggered/I");
  MuTrigMC->Branch("l2DiMuIsoTriggered",&l2DiMuIsoTriggered,"l2DiMuIsoTriggered/I");
  MuTrigMC->Branch("l3DiMuIsoPreTriggered",&l3DiMuIsoPreTriggered,"l3DiMuIsoPreTriggered/I");
  MuTrigMC->Branch("l3DiMuIsoTriggered",&l3DiMuIsoTriggered,"l3DiMuIsoTriggered/I"); 
  MuTrigMC->Branch("triggerDecisions",&triggerDecisions);
  MuTrigMC->Branch("triggerNames",&triggerNames);
  */
  // number of muons at each level
  MuTrigMC->Branch("nMu",&nMu,"nMu/I");
  MuTrigMC->Branch("nL2",&nL2,"nL2/I");
  MuTrigMC->Branch("nL3",&nL3,"nL3/I");
  MuTrigMC->Branch("nTkTracks",&nTkTracks,"nTkTracks/I");
  MuTrigMC->Branch("nL3Cands",&nL3Cands,"nL3Cands/I");
  MuTrigMC->Branch("muAllGlobalMuons",&muAllGlobalMuons);
  MuTrigMC->Branch("muAllStandAloneMuons",&muAllStandAloneMuons);
  MuTrigMC->Branch("muAllTrackerMuons",&muAllTrackerMuons);
  MuTrigMC->Branch("muTrackerMuonArbitrated",&muTrackerMuonArbitrated);
  MuTrigMC->Branch("muAllArbitrated",&muAllArbitrated);
  MuTrigMC->Branch("muGlobalMuonPromptTight",&muGlobalMuonPromptTight);
  MuTrigMC->Branch("muTMLastStationLoose",&muTMLastStationLoose);
  MuTrigMC->Branch("muTMLastStationTight",&muTMLastStationTight);
  MuTrigMC->Branch("muTM2DCompatibilityLoose",&muTM2DCompatibilityLoose);
  MuTrigMC->Branch("muTM2DCompatibilityTight",&muTM2DCompatibilityTight);
  MuTrigMC->Branch("muTMOneStationLoose",&muTMOneStationLoose);
  MuTrigMC->Branch("muTMOneStationTight",&muTMOneStationTight);
  MuTrigMC->Branch("muTMLastStationOptimizedLowPtLoose",&muTMLastStationOptimizedLowPtLoose);
  MuTrigMC->Branch("muTMLastStationOptimizedLowPtTight",&muTMLastStationOptimizedLowPtTight);
  MuTrigMC->Branch("muGMTkChiCompatibility",&muGMTkChiCompatibility);
  MuTrigMC->Branch("muGMStaChiCompatibility",&muGMStaChiCompatibility);
  MuTrigMC->Branch("muGMTkKinkTight",&muGMTkKinkTight);
  MuTrigMC->Branch("muCaloCompatibility",&muCaloCompatibility);
  MuTrigMC->Branch("muSegmentCompatibility",&muSegmentCompatibility);
  MuTrigMC->Branch("muTrkKink",&muTrkKink);
  MuTrigMC->Branch("muGlbKink",&muGlbKink);
  MuTrigMC->Branch("muTrkRelChi2",&muTrkRelChi2);
  MuTrigMC->Branch("muStaRelChi2",&muStaRelChi2);
  //  MuTrigMC->Branch("muIso03Valid",&muIso03Valid);
  //  MuTrigMC->Branch("muIso03sumPt",&muIso03sumPt);
  //  MuTrigMC->Branch("muIso03emEt",&muIso03emEt);
  //  MuTrigMC->Branch("muIso03hadEt",&muIso03hadEt);
  //  MuTrigMC->Branch("muIso03nTracks",&muIso03nTracks);
  //  MuTrigMC->Branch("muIso03trackerVetoPt",&muIso03trackerVetoPt);
  MuTrigMC->Branch("muNumberOfChambers",&muNumberOfChambers);
  MuTrigMC->Branch("muNumberOfMatches",&muNumberOfMatches);
  MuTrigMC->Branch("muStationMask",&muStationMask);
  MuTrigMC->Branch("muNCSCSeg",&muNCSCSeg);
  MuTrigMC->Branch("muNDTSeg",&muNDTSeg);
  MuTrigMC->Branch("muNRPCSeg",&muNRPCSeg);
  //MuTrigMC->Branch("nL3Seeds",&nL3Seeds,"nL3Seeds/I");
  // L3 muon information: the basics
  MuTrigMC->Branch("l3P",&l3P);
  MuTrigMC->Branch("l3Px",&l3Px);
  MuTrigMC->Branch("l3Py",&l3Py);
  MuTrigMC->Branch("l3Pz",&l3Pz);
  MuTrigMC->Branch("l3Pt",&l3Pt);
  MuTrigMC->Branch("l3PtError",&l3PtError);
  MuTrigMC->Branch("l3Pt90",&l3Pt90);
  MuTrigMC->Branch("l3Eta",&l3Eta);
  MuTrigMC->Branch("l3EtaError",&l3EtaError);
  MuTrigMC->Branch("l3Phi",&l3Phi);
  MuTrigMC->Branch("l3PhiError",&l3PhiError);
  MuTrigMC->Branch("l3D0",&l3D0);
  MuTrigMC->Branch("l3D0Error",&l3D0Error);
  MuTrigMC->Branch("l3NHits",&l3NHits);
  MuTrigMC->Branch("l3Charge",&l3Charge);
  MuTrigMC->Branch("l3Chi2",&l3Chi2);
  MuTrigMC->Branch("l3Ndof",&l3Ndof);
  // L3 Muon Isolation quantities
  MuTrigMC->Branch("l3CalIsoDeposit",&l3CalIsoDeposit);
  MuTrigMC->Branch("l3TrackIsoDeposit",&l3TrackIsoDeposit);
  MuTrigMC->Branch("l3IsoTrackDR",&l3IsoTrackDR);
  MuTrigMC->Branch("l3IsoTrackDRMinDelPt",&l3IsoTrackDRMinDelPt);
  // L3 Muon Track fitting parameters (with phi already declared above)
  MuTrigMC->Branch("l3Dsz",&l3Dsz);
  MuTrigMC->Branch("l3DszError",&l3DszError);
  MuTrigMC->Branch("l3Dxy",&l3Dxy);
  MuTrigMC->Branch("l3DxyError",&l3DxyError);
  MuTrigMC->Branch("l3Lambda",&l3Lambda);
  MuTrigMC->Branch("l3LambdaError",&l3LambdaError);
  MuTrigMC->Branch("l3Qoverp",&l3Qoverp);
  MuTrigMC->Branch("l3QoverpError",&l3QoverpError);
  MuTrigMC->Branch("l3ErrorMatrix",&l3ErrorMatrix);
  MuTrigMC->Branch("l3DetIds",&l3DetIds);
  MuTrigMC->Branch("l3SubdetIds",&l3SubdetIds);
  MuTrigMC->Branch("l3Component",&l3Component);
  MuTrigMC->Branch("l3NMuHits",&l3NMuHits);
  MuTrigMC->Branch("l3MuStationNumber",&l3MuStationNumber);
  MuTrigMC->Branch("l3RecHitsStatus",&l3RecHitsStatus);
  MuTrigMC->Branch("l3RecHitsX",&l3RecHitsX);
  MuTrigMC->Branch("l3RecHitsY",&l3RecHitsY);
  MuTrigMC->Branch("l3RecHitsZ",&l3RecHitsZ);
  MuTrigMC->Branch("l3RecHitsXTM",&l3RecHitsXTM);
  MuTrigMC->Branch("l3RecHitsYTM",&l3RecHitsYTM);
  MuTrigMC->Branch("l3RecHitsZTM",&l3RecHitsZTM);
  MuTrigMC->Branch("l3RecHitsXTSOS",&l3RecHitsXTSOS);
  MuTrigMC->Branch("l3RecHitsYTSOS",&l3RecHitsYTSOS);
  MuTrigMC->Branch("l3RecHitsZTSOS",&l3RecHitsZTSOS);
  MuTrigMC->Branch("l3RecHitsPhiTM",&l3RecHitsPhiTM);
  MuTrigMC->Branch("l3RecHitsErrorTM",&l3RecHitsErrorTM);
  MuTrigMC->Branch("l3RecHitsPhiTSOS",&l3RecHitsPhiTSOS);
  // Indices for L3<->L2
  MuTrigMC->Branch("indexL2SeedingL3",&indexL2SeedingL3);
  MuTrigMC->Branch("indexL3SeededFromL2",&indexL3SeededFromL2);
  MuTrigMC->Branch("l2SeedsL3",&l2SeedsL3);

  // The muon error matrix
  MuTrigMC->Branch("muonErrorMatrix",&muonErrorMatrix);

  // information from hltL3TkTracksFromL2
  MuTrigMC->Branch("tkTrackP",&tkTrackP);
  MuTrigMC->Branch("tkTrackPx",&tkTrackPx);
  MuTrigMC->Branch("tkTrackPy",&tkTrackPy);
  MuTrigMC->Branch("tkTrackPz",&tkTrackPz);
  MuTrigMC->Branch("tkTrackPt",&tkTrackPt);
  MuTrigMC->Branch("tkTrackPtError",&tkTrackPtError);
  MuTrigMC->Branch("tkTrackEta",&tkTrackEta);
  MuTrigMC->Branch("tkTrackEtaError",&tkTrackEtaError);
  MuTrigMC->Branch("tkTrackPhi",&tkTrackPhi);
  MuTrigMC->Branch("tkTrackPhiError",&tkTrackPhiError);
  MuTrigMC->Branch("tkTrackD0",&tkTrackD0);
  MuTrigMC->Branch("tkTrackD0Error",&tkTrackD0Error);
  MuTrigMC->Branch("tkTrackNHits",&tkTrackNHits);
  MuTrigMC->Branch("tkTrackCharge",&tkTrackCharge);
  MuTrigMC->Branch("tkTrackChi2",&tkTrackChi2);
  MuTrigMC->Branch("tkTrackNdof",&tkTrackNdof);
  // L3TRACK Muon Track fitting parameters (with phi already declared above)
  MuTrigMC->Branch("tkTrackDsz",&tkTrackDsz);
  MuTrigMC->Branch("tkTrackDszError",&tkTrackDszError);
  MuTrigMC->Branch("tkTrackDxy",&tkTrackDxy);
  MuTrigMC->Branch("tkTrackDxyError",&tkTrackDxyError);
  MuTrigMC->Branch("tkTrackLambda",&tkTrackLambda);
  MuTrigMC->Branch("tkTrackLambdaError",&tkTrackLambdaError);
  MuTrigMC->Branch("tkTrackQoverp",&tkTrackQoverp);
  MuTrigMC->Branch("tkTrackQoverpError",&tkTrackQoverpError);
  MuTrigMC->Branch("tkTrackErrorMatrix",&tkTrackErrorMatrix);
  MuTrigMC->Branch("tkTrackDetIds",&tkTrackDetIds);
  MuTrigMC->Branch("tkTrackSubdetIds",&tkTrackSubdetIds);
  MuTrigMC->Branch("tkTrackRecHitsStatus",&tkTrackRecHitsStatus);
  MuTrigMC->Branch("tkTrackRecHitsX",&tkTrackRecHitsX);
  MuTrigMC->Branch("tkTrackRecHitsY",&tkTrackRecHitsY);
  MuTrigMC->Branch("tkTrackRecHitsZ",&tkTrackRecHitsZ);

  // L2 muon information: the basics
  MuTrigMC->Branch("l2P",&l2P);
  MuTrigMC->Branch("l2Px",&l2Px);
  MuTrigMC->Branch("l2Py",&l2Py);
  MuTrigMC->Branch("l2Pz",&l2Pz);
  MuTrigMC->Branch("l2Pt",&l2Pt);
  MuTrigMC->Branch("l2PtError",&l2PtError);
  MuTrigMC->Branch("l2Pt90",&l2Pt90);
  MuTrigMC->Branch("l2Eta",&l2Eta);
  MuTrigMC->Branch("l2EtaError",&l2EtaError);
  MuTrigMC->Branch("l2Phi",&l2Phi);
  MuTrigMC->Branch("l2PhiError",&l2PhiError);
  MuTrigMC->Branch("l2D0",&l2D0);
  MuTrigMC->Branch("l2D0Error",&l2D0Error);
  MuTrigMC->Branch("l2NHits",&l2NHits);
  MuTrigMC->Branch("l2Charge",&l2Charge);
  MuTrigMC->Branch("l2Chi2",&l2Chi2);
  MuTrigMC->Branch("l2Ndof",&l2Ndof);
  MuTrigMC->Branch("l2DetIds",&l2DetIds);
  MuTrigMC->Branch("l2SubdetIds",&l2SubdetIds);
  MuTrigMC->Branch("l2Component",&l2Component);
  MuTrigMC->Branch("l2NMuHits",&l2NMuHits);
  MuTrigMC->Branch("l2MuStationNumber",&l2MuStationNumber);
  MuTrigMC->Branch("l2RecHitsStatus",&l2RecHitsStatus);
  MuTrigMC->Branch("l2RecHitsX",&l2RecHitsX);
  MuTrigMC->Branch("l2RecHitsY",&l2RecHitsY);
  MuTrigMC->Branch("l2RecHitsZ",&l2RecHitsZ);
  // L2 Muon Isolation quantities
  MuTrigMC->Branch("l2CalIsoDeposit",&l2CalIsoDeposit);
  // L2 Muon Track fitting parameters (with phi already declared above)
  MuTrigMC->Branch("l2Dsz",&l2Dsz);
  MuTrigMC->Branch("l2DszError",&l2DszError);
  MuTrigMC->Branch("l2Dxy",&l2Dxy);
  MuTrigMC->Branch("l2DxyError",&l2DxyError);
  MuTrigMC->Branch("l2Lambda",&l2Lambda);
  MuTrigMC->Branch("l2LambdaError",&l2LambdaError);
  MuTrigMC->Branch("l2Qoverp",&l2Qoverp);
  MuTrigMC->Branch("l2QoverpError",&l2QoverpError);
  MuTrigMC->Branch("l2ErrorMatrix",&l2ErrorMatrix);

  // Specific to MuTrigMC
  MuTrigMC->Branch("l3IsAssociated",&l3IsAssociated);
  MuTrigMC->Branch("l3ParentID",&l3ParentID);
  MuTrigMC->Branch("l3MotherBinNumber",&l3MotherBinNumber);
  MuTrigMC->Branch("l3AssociationVar",&l3AssociationVar);
  MuTrigMC->Branch("l3AssociationPdgId",&l3AssociationPdgId);
  MuTrigMC->Branch("l3AssociationMyBit",&l3AssociationMyBit);
  MuTrigMC->Branch("l3AssociationVtxX",&l3AssociationVtxX);
  MuTrigMC->Branch("l3AssociationVtxY",&l3AssociationVtxY);
  MuTrigMC->Branch("l3AssociationVtxZ",&l3AssociationVtxZ);
  MuTrigMC->Branch("l3AssociatedSimMuonIndex",&l3AssociatedSimMuonIndex);
  MuTrigMC->Branch("l3AssociatedSimMuonPt",&l3AssociatedSimMuonPt);
  MuTrigMC->Branch("l3AssociatedSimMuonEta",&l3AssociatedSimMuonEta);
  MuTrigMC->Branch("l3AssociatedSimMuonPhi",&l3AssociatedSimMuonPhi);
  MuTrigMC->Branch("l3AssociatedSimMuonNHits",&l3AssociatedSimMuonNHits);
  MuTrigMC->Branch("l3AssociatedSimMuonDetIds",&l3AssociatedSimMuonDetIds);
  MuTrigMC->Branch("l3AssociatedSimMuonNMuHits",&l3AssociatedSimMuonNMuHits);
  MuTrigMC->Branch("l3AssociatedSimMuonMuStationNumber",&l3AssociatedSimMuonMuStationNumber);

  MuTrigMC->Branch("l3AssociatedSimMuonDsz",&l3AssociatedSimMuonDsz);
  MuTrigMC->Branch("l3AssociatedSimMuonDxy",&l3AssociatedSimMuonDxy);
  MuTrigMC->Branch("l3AssociatedSimMuonLambda",&l3AssociatedSimMuonLambda);
  MuTrigMC->Branch("l3AssociatedSimMuonQoverP",&l3AssociatedSimMuonQoverP);

  MuTrigMC->Branch("tkTrackIsAssociated",&tkTrackIsAssociated);
  MuTrigMC->Branch("tkTrackParentID",&tkTrackParentID);
  MuTrigMC->Branch("tkTrackMotherBinNumber",&tkTrackMotherBinNumber);
  MuTrigMC->Branch("tkTrackAssociationVar",&tkTrackAssociationVar);
  MuTrigMC->Branch("tkTrackAssociationPdgId",&tkTrackAssociationPdgId);
  MuTrigMC->Branch("tkTrackAssociationMyBit",&tkTrackAssociationMyBit);
  MuTrigMC->Branch("tkTrackAssociationVtxX",&tkTrackAssociationVtxX);
  MuTrigMC->Branch("tkTrackAssociationVtxY",&tkTrackAssociationVtxY);
  MuTrigMC->Branch("tkTrackAssociationVtxZ",&tkTrackAssociationVtxZ);
  MuTrigMC->Branch("tkTrackAssociatedSimMuonPt",&tkTrackAssociatedSimMuonPt);
  MuTrigMC->Branch("tkTrackAssociatedSimMuonEta",&tkTrackAssociatedSimMuonEta);
  MuTrigMC->Branch("tkTrackAssociatedSimMuonPhi",&tkTrackAssociatedSimMuonPhi);
  MuTrigMC->Branch("tkTrackAssociatedSimMuonNHits",&tkTrackAssociatedSimMuonNHits);
  MuTrigMC->Branch("tkTrackAssociatedSimMuonDetIds",&tkTrackAssociatedSimMuonDetIds);

  MuTrigMC->Branch("tkTrackAssociatedSimMuonDsz",&tkTrackAssociatedSimMuonDsz);
  MuTrigMC->Branch("tkTrackAssociatedSimMuonDxy",&tkTrackAssociatedSimMuonDxy);
  MuTrigMC->Branch("tkTrackAssociatedSimMuonLambda",&tkTrackAssociatedSimMuonLambda);
  MuTrigMC->Branch("tkTrackAssociatedSimMuonQoverP",&tkTrackAssociatedSimMuonQoverP);

  MuTrigMC->Branch("l2IsAssociated",&l2IsAssociated);
  MuTrigMC->Branch("l2ParentID",&l2ParentID);
  MuTrigMC->Branch("l2MotherBinNumber",&l2MotherBinNumber);
  MuTrigMC->Branch("l2AssociationVar",&l2AssociationVar);
  MuTrigMC->Branch("l2AssociationPdgId",&l2AssociationPdgId);
  MuTrigMC->Branch("l2AssociationMyBit",&l2AssociationMyBit);
  MuTrigMC->Branch("l2AssociationVtxX",&l2AssociationVtxX);
  MuTrigMC->Branch("l2AssociationVtxY",&l2AssociationVtxY);
  MuTrigMC->Branch("l2AssociationVtxZ",&l2AssociationVtxZ);
  MuTrigMC->Branch("l2AssociatedSimMuonPt",&l2AssociatedSimMuonPt);
  MuTrigMC->Branch("l2AssociatedSimMuonEta",&l2AssociatedSimMuonEta);
  MuTrigMC->Branch("l2AssociatedSimMuonPhi",&l2AssociatedSimMuonPhi);
  MuTrigMC->Branch("l2AssociatedSimMuonNHits",&l2AssociatedSimMuonNHits);
  MuTrigMC->Branch("l2AssociatedSimMuonDetIds",&l2AssociatedSimMuonDetIds);
  MuTrigMC->Branch("l2AssociatedSimMuonNMuHits",&l2AssociatedSimMuonNMuHits);
  MuTrigMC->Branch("l2AssociatedSimMuonMuStationNumber",&l2AssociatedSimMuonMuStationNumber);

  MuTrigMC->Branch("l2AssociatedSimMuonDsz",&l2AssociatedSimMuonDsz);
  MuTrigMC->Branch("l2AssociatedSimMuonDxy",&l2AssociatedSimMuonDxy);
  MuTrigMC->Branch("l2AssociatedSimMuonLambda",&l2AssociatedSimMuonLambda);
  MuTrigMC->Branch("l2AssociatedSimMuonQoverP",&l2AssociatedSimMuonQoverP);

  MuTrigMC->Branch("nSimMuon",&nSimMuon,"nSimMuon/I");
  MuTrigMC->Branch("simMuonParentID",&simMuonParentID);
  MuTrigMC->Branch("simMuonMotherBinNumber",&simMuonMotherBinNumber);
  MuTrigMC->Branch("simMuonPt",&simMuonPt);
  MuTrigMC->Branch("simMuonEta",&simMuonEta);
  MuTrigMC->Branch("simMuonPhi",&simMuonPhi);
  MuTrigMC->Branch("simMuonMyBit",&simMuonMyBit);
  MuTrigMC->Branch("simMuonVtxX",&simMuonVtxX);
  MuTrigMC->Branch("simMuonVtxY",&simMuonVtxY);
  MuTrigMC->Branch("simMuonVtxZ",&simMuonVtxZ);
  MuTrigMC->Branch("simMuonNHits",&simMuonNHits);
  MuTrigMC->Branch("simMuonDetIds",&simMuonDetIds);
  MuTrigMC->Branch("simMuonNMuHits",&simMuonNMuHits);
  MuTrigMC->Branch("simMuonMuStationNumber",&simMuonMuStationNumber);

  MuTrigMC->Branch("simMuonDsz",&simMuonDsz);
  MuTrigMC->Branch("simMuonDxy",&simMuonDxy);
  MuTrigMC->Branch("simMuonLambda",&simMuonLambda);
  MuTrigMC->Branch("simMuonQoverP",&simMuonQoverP);

  MuTrigMC->Branch("simToL3Associated",&simToL3Associated);
  MuTrigMC->Branch("simToL3AssociationVar",&simToL3AssociationVar);
  MuTrigMC->Branch("simToL3RecoIndex",&simToL3RecoIndex);
  MuTrigMC->Branch("simToTkAssociated",&simToTkAssociated);
  MuTrigMC->Branch("simToTkAssociationVar",&simToTkAssociationVar);
  MuTrigMC->Branch("simToTkRecoIndex",&simToTkRecoIndex);
  MuTrigMC->Branch("simToL2Associated",&simToL2Associated);
  MuTrigMC->Branch("simToL2AssociationVar",&simToL2AssociationVar);
  MuTrigMC->Branch("simToL2RecoIndex",&simToL2RecoIndex);

  edm::LogInfo("MuonRecoTreeUtility")<<"beginJob executed.  Problems not from here.";

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonRecoTreeUtility::endJob() {

  edm::LogInfo("MuonRecoTreeUtility")<<"Starting to write the trees.  Changing to directory";
  theFile->cd();
  edm::LogInfo("MuonRecoTreeUtility")<<"Starting to write the trees.  Writing MuTrigData";
  MuTrigData->Write();
  edm::LogInfo("MuonRecoTreeUtility")<<"Starting to write the trees.  Writing MuTrigMC";
  MuTrigMC->Write();
  edm::LogInfo("MuonRecoTreeUtility")<<"Finished writing trees.  Closing file.";
  theFile->Close();
  edm::LogInfo("MuonRecoTreeUtility")<<"All done.  Nothing left to do.";

}

unsigned int MuonRecoTreeUtility::getBit(const TrackingParticleRef& simRef) const
{
  /*  
  LogDebug("SpecialBit") << "\nRecoTree\nTrackingParticle " << simRef->pdgId()
			 << " with pt= " << sqrt(simRef->momentum().perp2())
			 << " at r " << sqrt(simRef->vertex().x()*simRef->vertex().x() + simRef->vertex().y()*simRef->vertex().y()) << " and z " << simRef->vertex().z() << "\n";
  */

  unsigned int thisBit = noBit;
  
  if ( tpSelector_primary(*simRef) ) {
    thisBit = primaryMuon;
    LogDebug("SpecialBit") << "b1";
  }
  else if ( tpSelector_silicon(*simRef) && !(isPrimaryMuon(thisBit)) ){
    thisBit = siliconMuon;
    LogDebug("SpecialBit") << "b2";
  }
  else if ( tpSelector_calConversion(*simRef) && ! (isSiliconMuon(thisBit) ||  isPrimaryMuon(thisBit)) ) {
    thisBit = calConversionMuon;
    LogDebug("SpecialBit") << "b3";
  }
  else if ( fabs(simRef->pdgId())==13 && simRef->pt() >= 0.9 && ! (isPrimaryMuon(thisBit) ||  isSiliconMuon(thisBit) ||  isCalConversionMuon(thisBit)) ) {
    thisBit = otherMuon;
    LogDebug("SpecialBit") << "b4";
  }
  
  return thisBit; 
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonRecoTreeUtility);
