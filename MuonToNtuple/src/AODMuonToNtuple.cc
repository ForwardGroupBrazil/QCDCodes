/** \class AODMuonToNtuple
 *  $Date: 2009/07/25 16:47:13 $
 *  $Revision: 1.3 $
 *  \author Chang Liu   -  Purdue University <Chang.Liu@cern.ch>
 *
 *  1st modified by Hwidong Yoo (hdyoo@cern.ch)
 *  Sep. 15, 2008
 *
 *  version 0.1 by Hwidong Yoo (hdyoo@cern.ch)
 *  Oct. 20, 2008
 */


#include "UserCode/MuonToNtuple/interface/AODMuonToNtuple.h"
// system include files
#include <memory>


//#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/Vertexing.h"
//#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// more info
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

// trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerRecord.h"

#include "RecoMuon/GlobalTrackingTools/interface/GlobalMuonRefitter.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/PatternTools/interface/MeasurementEstimator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackHistory/interface/TrackClassifier.h"
#include "SimTracker/TrackHistory/interface/TrackCategories.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include <TrackingTools/PatternTools/interface/TrajectoryMeasurement.h>
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
#include <TrackingTools/PatternTools/interface/Trajectory.h>
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include <DataFormats/TrackingRecHit/interface/TrackingRecHit.h>
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "DataFormats/DetId/interface/DetId.h"

#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFrame.h>
#include <TMath.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TVector3.h>

//
// class decleration
//

using namespace std;
using namespace reco;
using namespace edm;


//
AODMuonToNtuple::AODMuonToNtuple(const edm::ParameterSet& iConfig)
{
  nEvt = 0;
  LogDebug("NTUPLE");
  theRootFileName = iConfig.getParameter<string>("out");
  theMuonLabel = iConfig.getParameter<edm::InputTag>("Muon");
  theMETLabel = iConfig.getParameter<edm::InputTag>("MET");
  theJetLabel = iConfig.getParameter<edm::InputTag>("Jet");
  isMC = iConfig.getParameter<bool>("isMC");
  theCrossSection = iConfig.getParameter<double>("CrossSection");
  theFilterEfficiency = iConfig.getParameter<double>("FilterEfficiency");
  theTotalNevents = iConfig.getParameter<double>("TotalNevents");
  theBDiscriminant = iConfig.getParameter<double>("BDiscriminant");
  LogDebug("NTUPLE");
  //aaa classifier_ = new TrackClassifier(iConfig);
  // service parameters
  ParameterSet serviceParameters = iConfig.getParameter<ParameterSet>("ServiceParameters");

  // the services
  theService = new MuonServiceProxy(serviceParameters);
    LogDebug("NTUPLE");
  // TrackRefitter parameters
  ParameterSet refitterParameters = iConfig.getParameter<ParameterSet>("RefitterParameters");
  theRefitter = new GlobalMuonRefitter(refitterParameters, theService);

  double maxChi2 = iConfig.getParameter<double>("MaxChi2");
  double nSigma = iConfig.getParameter<double>("nSigma");
  LogDebug("NTUPLE");
  theEstimator = new Chi2MeasurementEstimator(maxChi2,nSigma);
  LogDebug("NTUPLE");
  /*
  ParameterSet trackTransformerPSet = iConfig.getParameter<ParameterSet>("TrackTransformer");
  //trackTransformerPSet.addParameter<string>("Propagator",TransformerOutPropagator);
  theTrackTransformer = new TrackTransformer(trackTransformerPSet);
  */
}


AODMuonToNtuple::~AODMuonToNtuple()
{
  if (theService) delete theService;
  if (theRefitter) delete theRefitter;
  if (theEstimator) delete theEstimator;
  //aaa if (classifier_) delete classifier_;

}


//
// member functions
//

// ------------ method called to for each event  ------------
void AODMuonToNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  LogDebug("NTUPLE");
  theService->update(iSetup);
  
  theRefitter->setEvent(iEvent);
  
  theRefitter->setServices(theService->eventSetup());

  //aaa classifier_->newEvent(iEvent, iSetup);
  LogDebug("NTUPLE");

   // initialize for ntuple variables
   nPair = -1; GENnPair = -1; 
   weight = -1; MET = -1;
   Nmuons = Njets = Nbtagged = NbtaggedCloseMuon = -1;
   PVtrackSize = -1;
   PVchi2 = -1;
   PVndof = -1;
   PVnormalizedChi2 = -1;
   PVx = -1000;
   PVy = -1000;
   PVz = -1000;

   // trigger lists
   _HLT_L1Mu = -1;
   _HLT_L1MuOpen = -1;
   _HLT_L2Mu9 = -1;
   _HLT_IsoMu9 = -1;   _HLT_IsoMu11 = -1;
   _HLT_IsoMu13 = -1;
   _HLT_IsoMu15 = -1;
   _HLT_Mu3 = -1;
   _HLT_Mu5 = -1;
   _HLT_Mu7 = -1;
   _HLT_Mu9 = -1;
   _HLT_Mu11 = -1;
   _HLT_Mu13 = -1;
   _HLT_Mu15 = -1;
   _HLT_Mu15_L1Mu7 = -1;
   _HLT_Mu15_Vtx2cm = -1;
   _HLT_Mu15_Vtx2mm = -1;
   _HLT_DoubleIsoMu3 = -1;
   _HLT_DoubleMu3 = -1;
   _HLT_DoubleMu3_Vtx2cm = -1;
   _HLT_DoubleMu3_Vtx2mm = -1;
   _HLT_DoubleMu3_JPsi = -1;
   _HLT_DoubleMu3_Upsilon = -1;
   _HLT_DoubleMu7_Z = -1;
   _HLT_DoubleMu3_SameSign = -1;
   _HLT_DoubleMu3_Psi2S = -1;

   for( int i = 0; i < MPSIZE; i++ ) {

       // Jet
       JETbDiscriminant[i] = JETcharge[i] = JETpt[i] = JETeta[i] = JETphi[i] = -100;
       JETflavour[i] = JETntracks[i] = -100;

       // muon type
       Muon1_muonType[i] = -1;

       // trigger
       Muon1_nTrig[i] = Muon1_triggerObjectType[i] = Muon1_filterName[i] = -1;

       // Muon kinematics
       InvMass[i] = Angle[i] = isOppSign[i] = MuonPairIndex[i] = isHMassPair[i] = -1;
       DeltaR[i] = Dalpha[i] = -100;
       vtxChi2[i] = -100;
       vtxPositionX[i] = vtxPositionY[i] = vtxPositionZ[i] = -99999;
       vtxXerror[i] = vtxYerror[i] = vtxZerror[i] = -99999;

       Muon1_phi[i] = Muon1_eta[i] = Muon1_pT[i] = -100;
       Muon1_Px[i] = Muon1_Py[i] = Muon1_Pz[i] = -100;
       Muon1_Glb_phi[i] = Muon1_Glb_eta[i] = Muon1_Glb_pT[i] = -100;
       Muon1_Glb_Px[i] = Muon1_Glb_Py[i] = Muon1_Glb_Pz[i] = -100;
       Muon1_Sta_phi[i] = Muon1_Sta_eta[i] = Muon1_Sta_pT[i] = -100;
       Muon1_Sta_Px[i] = Muon1_Sta_Py[i] = Muon1_Sta_Pz[i] = -100;
       Muon1_Trk_phi[i] = Muon1_Trk_eta[i] = Muon1_Trk_pT[i] = -100;
       Muon1_Trk_Px[i] = Muon1_Trk_Py[i] = Muon1_Trk_Pz[i] = -100;
       Muon1_trkiso[i] = Muon1_hcaliso[i] = Muon1_ecaliso[i] = Muon1_chi2dof[i] = -100;
       Muon1_nChambers[i] = Muon1_nMatches[i] = Muon1_stationMask[i] = Muon1_nSegments[i] =  -1;

       Muon1_charge[i] = Muon1_nhits[i] = -100;
       Muon1_Glb_charge[i] = Muon1_Glb_nhits[i] = -100;
       Muon1_Sta_charge[i] = Muon1_Sta_nhits[i] = -100;
       Muon1_Trk_charge[i] = Muon1_Trk_nhits[i] = -100;

       Muon1_qoverp[i] = Muon1_theta[i] = Muon1_lambda[i] = -100;
       Muon1_dxy[i] = Muon1_d0[i] = Muon1_dsz[i] = Muon1_dz[i] = -100;
       Muon1_vx[i] = Muon1_vy[i] = Muon1_vz[i] = -100;
       Muon1_dxyBS[i] = Muon1_dszBS[i] = Muon1_dzBS[i] = -100;

       Muon1_Glb_qoverp[i] = Muon1_Glb_theta[i] = Muon1_Glb_lambda[i] = -100;
       Muon1_Glb_dxy[i] = Muon1_Glb_d0[i] = Muon1_Glb_dsz[i] = Muon1_Glb_dz[i] = -100;
       Muon1_Glb_vx[i] = Muon1_Glb_vy[i] = Muon1_Glb_vz[i] = -100;
       Muon1_Glb_trkKink[i] = Muon1_Glb_glbKink[i] = -100;
       Muon1_Glb_rChi2Sta[i] = Muon1_Glb_rChi2Trk[i] = -100;
       Muon1_Glb_dxyBS[i] = Muon1_Glb_dszBS[i] = Muon1_Glb_dzBS[i] = -100;

       Muon1_Sta_qoverp[i] = Muon1_Sta_theta[i] = Muon1_Sta_lambda[i] = -100;
       Muon1_Sta_dxy[i] = Muon1_Sta_d0[i] = Muon1_Sta_dsz[i] = Muon1_Sta_dz[i] = -100;
       Muon1_Sta_vx[i] = Muon1_Sta_vy[i] = Muon1_Sta_vz[i] = -100;
       Muon1_Sta_dxyBS[i] = Muon1_Sta_dszBS[i] = Muon1_Sta_dzBS[i] = -100;

       Muon1_Trk_qoverp[i] = Muon1_Trk_theta[i] = Muon1_Trk_lambda[i] = -100;
       Muon1_Trk_dxy[i] = Muon1_Trk_d0[i] = Muon1_Trk_dsz[i] = Muon1_Trk_dz[i] = -100;
       Muon1_Trk_vx[i] = Muon1_Trk_vy[i] = Muon1_Trk_vz[i] = -100;
       Muon1_Trk_dxyBS[i] = Muon1_Trk_dszBS[i] = Muon1_Trk_dzBS[i] = -100;

       //       Muon2_qoverp[i] = Muon2_theta[i] = Muon2_lambda[i] = -100;
       //       Muon2_dxy[i] = Muon2_d0[i] = Muon2_dsz[i] = Muon2_dz[i] = -100;
       //       Muon2_vx[i] = Muon2_vy[i] = Muon2_vz[i] = -100;
       //       Muon2_dxyBS[i] = Muon2_dszBS[i] = Muon2_dzBS[i] = -100;

       Muon1_MCtruth_pT[i] = Muon1_MCtruth_eta[i] = Muon1_MCtruth_phi[i] = -100;
       Muon1_MCtruth_Px[i] = Muon1_MCtruth_Py[i] = Muon1_MCtruth_Pz[i] = -100;
       Muon1_MCtruth_charge[i] = Muon1_MCtruth_mother[i] =-100;

       // GEN Muon
       GENInvMass[i] = GENAngle[i] = GENisOppSign[i] = -1;
       GENMuon1_phi[i] = GENMuon1_eta[i] = GENMuon1_pT[i] = GENMuon1_mother[i] = -100;
       GENMuon2_phi[i] = GENMuon2_eta[i] = GENMuon2_pT[i] = GENMuon2_mother[i] = -100;
       GENMuon1_Px[i] = GENMuon1_Py[i] = GENMuon1_Pz[i] = -100;
       GENMuon2_Px[i] = GENMuon2_Py[i] = GENMuon2_Pz[i] = -100;
       GENMuon1_charge[i] = GENMuon2_charge[i] = -100;
   }

   const double _intLumi = 100; // set to 100 pb-1 integrated luminosity
   using namespace edm;

   LogDebug("NTUPLE") << "reading event " << iEvent.id() ;
   nEvt++;
   // run number & event number
   runNum = iEvent.id().run();
   evtNum = iEvent.id().event();

   // call AOD objects
   edm::Handle<edm::View<reco::Muon> > muonHandle;
   iEvent.getByLabel(theMuonLabel,muonHandle);

   edm::Handle<edm::View<pat::MET> > metHandle;
   iEvent.getByLabel(theMETLabel,metHandle);
   
   edm::Handle<edm::View<pat::Jet> > jetHandle;
   iEvent.getByLabel(theJetLabel,jetHandle);
   
   edm::Handle<edm::View<reco::Track> > trackHandle;
   iEvent.getByLabel("generalTracks", trackHandle);

   edm::Handle<edm::View<reco::Track> > standAloneHandle;
   iEvent.getByLabel("standAloneMuons:UpdatedAtVtx", standAloneHandle);

   edm::Handle<edm::View<reco::Track> > globalHandle;
   iEvent.getByLabel("globalMuons", globalHandle);

   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
   reco::BeamSpot beamSpot = (*beamSpotHandle);

   // Primary Vertex
   edm::Handle<reco::VertexCollection> pvHandle;
   iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
   const reco::VertexCollection vtx = *(pvHandle.product());

   if( vtx.size() > 2 ) LogDebug("NTUPLE") << "Reconstructed "<< vtx.size() << " vertices" ;
   if (vtx.size() > 0 ){
       PVtrackSize = vtx.front().tracksSize();
       PVchi2 = vtx.front().chi2();
       PVndof = vtx.front().ndof();
       PVnormalizedChi2 = vtx.front().normalizedChi2();
       PVx = vtx.front().x();
       PVy = vtx.front().y();
       PVz = vtx.front().z();
   }

   // only choose an event with more than two muons
   Nmuons = muonHandle->size();
   if (Nmuons < 2) return;

   
   //
   // trigger
   //
   /*
   const int nTrigName = 26;
   // put the whole muon HLT trigger lists defined in CMSSW_2_1_9
   string _trigs[nTrigName] = {"HLT_L1Mu", "HLT_L1MuOpen", "HLT_L2Mu9", 
   		"HLT_IsoMu9", "HLT_IsoMu11", "HLT_IsoMu13", "HLT_IsoMu15",
   		"HLT_Mu3", "HLT_Mu5", "HLT_Mu7", "HLT_Mu9", "HLT_Mu11", "HLT_Mu13", "HLT_Mu15", 
		"HLT_Mu15_L1Mu7", "HLT_Mu15_Vtx2cm", "HLT_Mu15_Vtx2mm", 
		"HLT_DoubleIsoMu3", "HLT_DoubleMu3", "HLT_DoubleMu3_Vtx2cm", "HLT_DoubleMu3_Vtx2mm", 
		"HLT_DoubleMu3_JPsi", "HLT_DoubleMu3_Upsilon", "HLT_DoubleMu7_Z",
		"HLT_DoubleMu3_SameSign", "HLT_DoubleMu3_Psi2S"};
   */
   const int nTrigName = 1;
   // put the whole muon HLT trigger lists defined in CMSSW_2_1_9
   string _trigs[nTrigName] = {"HLT_Mu9"} ;


/*
   const int nL1TrigName = 10;
   string _L1trigs[nTrigName] = {"L1_SingleMu0", "L1_SingleMu3", "L1_SingleMu5",
   		"L1_SingleMu7", "L1_SingleMu10", "L1_SingleMu14", "L1_SingleMu20",
       		"L1_SingleMu25", "L1_DoubleMu3", "L1_TripleMu3"};
		*/

   // L1 Trigger
   edm::Handle<L1GlobalTriggerReadoutRecord> gtReadoutRecord;
   iEvent.getByLabel(InputTag("gtDigis"), gtReadoutRecord);

   if( gtReadoutRecord.isValid() ) {
       // get Global Trigger finalOR and the decision word
       boost::uint16_t gtFinalOR = gtReadoutRecord->finalOR();
       //gtDecisionWordBeforeMask = gtReadoutRecord->decisionWord();
       //technicalTriggerWordBeforeMask = gtReadoutRecord->technicalTriggerWord();
       unsigned int m_numberDaqPartitions = 8;
       for( unsigned int iDaqPartition = 0; iDaqPartition < m_numberDaqPartitions; ++iDaqPartition ) {
	   bool gtDecision = static_cast<bool>(gtFinalOR & ( 1 << iDaqPartition ));
	   if( gtDecision ) {
	   }
       }
   }
   
   std::vector<std::string > MuonHLT;
   for( int i = 0; i < nTrigName; i++ ) {
       MuonHLT.push_back(_trigs[i]);
   }

   // read the whole HLT trigger lists fired in an event
   bool trigFired[nTrigName] = {false};
   Handle<TriggerResults> trigResult;
   iEvent.getByLabel(InputTag("TriggerResults","","HLT"), trigResult);
   int ntrigs = trigResult->size();   LogDebug("NTUPLE") << ntrigs;
   TriggerNames trigName;
   if( ntrigs != 0 ) {
       trigName.init(*trigResult);
       for( int itrig = 0; itrig != ntrigs; ++itrig) {
       //aaa2	   for( int itrigName = 0; itrigName < nTrigName; itrigName++ ) {
	   for( int itrigName = 0; itrigName < 0; itrigName++ ) {
	     // LogDebug("NTUPLE") << itrig << " " << itrigName;
	     if( trigResult->accept(trigName.triggerIndex(MuonHLT[itrigName])) ){ trigFired[itrigName] = true;}
	   }
       }
       _HLT_L1Mu = trigFired[0];
       _HLT_L1MuOpen = trigFired[1];
       _HLT_L2Mu9 = trigFired[2];
       _HLT_IsoMu9 = trigFired[3];
       _HLT_IsoMu11 = trigFired[4];
       _HLT_IsoMu13 = trigFired[5];
       _HLT_IsoMu15 = trigFired[6];
       _HLT_Mu3 = trigFired[7];
       _HLT_Mu5 = trigFired[8];
       _HLT_Mu7 = trigFired[9];
       _HLT_Mu9 = trigFired[10];
       _HLT_Mu11 = trigFired[11];
       _HLT_Mu13 = trigFired[12];
       _HLT_Mu15 = trigFired[13];
       _HLT_Mu15_L1Mu7 = trigFired[14];
       _HLT_Mu15_Vtx2cm = trigFired[15];
       _HLT_Mu15_Vtx2mm = trigFired[16];
       _HLT_DoubleIsoMu3 = trigFired[17];
       _HLT_DoubleMu3 = trigFired[18];
       _HLT_DoubleMu3_Vtx2cm = trigFired[19];
       _HLT_DoubleMu3_Vtx2mm = trigFired[20];
       _HLT_DoubleMu3_JPsi = trigFired[21];
       _HLT_DoubleMu3_Upsilon = trigFired[22];
       _HLT_DoubleMu7_Z = trigFired[23];
       _HLT_DoubleMu3_SameSign = trigFired[24];
       _HLT_DoubleMu3_Psi2S = trigFired[25];
   }

   //
   // select highest invariant mass pair
   //
   /* 
   double highest_invmass = -100;
   reco::Muon hMmuon1; reco::Muon hMmuon2;
   for(edm::View<reco::Muon>::const_iterator iMuon1 = muonHandle->begin(); iMuon1 != muonHandle->end(); ++iMuon1) {
       for(edm::View<reco::Muon>::const_iterator iMuon2 = muonHandle->begin(); iMuon2 != muonHandle->end(); ++iMuon2) {
	    if( iMuon1 <= iMuon2 ) continue;
	    if( !iMuon1->isGood() || !iMuon2->isGood() ) continue;

	    reco::NamedCompositeCandidate aDYcand;
    	    aDYcand.setP4(iMuon1->p4() + iMuon2->p4());
	    double invmass = aDYcand.mass();
	    if( invmass > highest_invmass ) {
		hMmuon1 = *iMuon1;
		hMmuon2 = *iMuon2;
	    }
       }
   }
   */

   //
   // MuonTree (loop over all muons)
   //
   int _nPair = 0;
   for(edm::View<reco::Muon>::const_iterator iMuon1 = muonHandle->begin(); iMuon1 != muonHandle->end(); ++iMuon1) {
     //aaa01 for(edm::View<reco::Muon>::const_iterator iMuon2 = muonHandle->begin(); iMuon2 != muonHandle->end(); ++iMuon2) {
     //if( iMuon1 <= iMuon2 ) continue;
	    // only select good muons
     //aaa2	    if( !iMuon1->isGood() ) continue;
	    if( !muon::isGoodMuon(*iMuon1,muon::AllArbitrated) ) continue;
	    LogDebug("NTUPLE");
	    // trigger matching in PAT
	    /*
	    Muon1_nTrig[_nPair] = iMuon1->triggerMatches().size();
	    if( Muon1_nTrig[_nPair] > 0 ) {
	    	// trigger object type as defined in DataFormats/HLTReco/interface/TriggerTypeDefs.h
	    	Muon1_triggerObjectType[_nPair] = iMuon1->triggerMatches()[0].triggerObjectType();
		string _filterName = iMuon1->triggerMatches()[0].filterName();
	    	if( _filterName == "hltMuLevel1PathL1Filtered" ) Muon1_filterName[_nPair] = 0;
	    	if( _filterName == "hltMuLevel1PathL1OpenFiltered" ) Muon1_filterName[_nPair] = 1;
	    	if( _filterName == "hltSingleMuLevel2NoIsoL2PreFiltered " ) Muon1_filterName[_nPair] = 2;
	    	if( _filterName == "hltSingleMuIsoL3IsoFiltered9" ) Muon1_filterName[_nPair] = 3;
	    	if( _filterName == "hltSingleMuIsoL3IsoFiltered" ) Muon1_filterName[_nPair] = 4;
	    	if( _filterName == "hltSingleMuIsoL3IsoFiltered13" ) Muon1_filterName[_nPair] = 5;
	    	if( _filterName == "hltSingleMuIsoL3IsoFiltered15" ) Muon1_filterName[_nPair] = 6;
	    	if( _filterName == "hltSingleMuPrescale3L3PreFiltered" ) Muon1_filterName[_nPair] = 7;
	    	if( _filterName == "hltSingleMuPrescale5L3PreFiltered" ) Muon1_filterName[_nPair] = 8;
	    	if( _filterName == "hltSingleMuPrescale77L3PreFiltered" ) Muon1_filterName[_nPair] = 9;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFiltered9" ) Muon1_filterName[_nPair] = 10;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFiltered11" ) Muon1_filterName[_nPair] = 11;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFiltered13" ) Muon1_filterName[_nPair] = 12;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFiltered15" ) Muon1_filterName[_nPair] = 13;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFiltered" ) Muon1_filterName[_nPair] = 14;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFilteredRelaxedVtx2cm" ) Muon1_filterName[_nPair] = 15;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFilteredRelaxedVtx2mm" ) Muon1_filterName[_nPair] = 16;
	    	if( _filterName == "hltDiMuonIsoL3IsoFiltered" ) Muon1_filterName[_nPair] = 17;
	    	if( _filterName == "hltDiMuonNoIsoL3PreFiltered" ) Muon1_filterName[_nPair] = 18;
	    	if( _filterName == "hltDiMuonNoIsoL3PreFilteredRelaxedVtx2cm" ) Muon1_filterName[_nPair] = 19;
	    	if( _filterName == "hltDiMuonNoIsoL3PreFilteredRelaxedVtx2mm" ) Muon1_filterName[_nPair] = 20;
	    	if( _filterName == "hltJpsiMML3Filtered" ) Muon1_filterName[_nPair] = 21;
	    	if( _filterName == "hltUpsilonMML3Filtered" ) Muon1_filterName[_nPair] = 22;
	    	if( _filterName == "hltZMML3Filtered" ) Muon1_filterName[_nPair] = 23;
	    	if( _filterName == "hltSameSignMuL3PreFiltered" ) Muon1_filterName[_nPair] = 24;
	    	if( _filterName == "hltPsi2SMML3Filtered" ) Muon1_filterName[_nPair] = 25;
	    }
	    */
	    /*
	    Muon2_nTrig[_nPair] = iMuon2->triggerMatches().size();
	    if( Muon2_nTrig[_nPair] > 0 ) {
	    	// trigger object type as defined in DataFormats/HLTReco/interface/TriggerTypeDefs.h
	    	Muon2_triggerObjectType[_nPair] = iMuon2->triggerMatches()[0].triggerObjectType();
		string _filterName = iMuon2->triggerMatches()[0].filterName();
	    	if( _filterName == "hltMuLevel1PathL1Filtered" ) Muon2_filterName[_nPair] = 0;
	    	if( _filterName == "hltMuLevel1PathL1OpenFiltered" ) Muon2_filterName[_nPair] = 1;
	    	if( _filterName == "hltSingleMuLevel2NoIsoL2PreFiltered " ) Muon2_filterName[_nPair] = 2;
	    	if( _filterName == "hltSingleMuIsoL3IsoFiltered9" ) Muon2_filterName[_nPair] = 3;
	    	if( _filterName == "hltSingleMuIsoL3IsoFiltered" ) Muon2_filterName[_nPair] = 4;
	    	if( _filterName == "hltSingleMuIsoL3IsoFiltered13" ) Muon2_filterName[_nPair] = 5;
	    	if( _filterName == "hltSingleMuIsoL3IsoFiltered15" ) Muon2_filterName[_nPair] = 6;
	    	if( _filterName == "hltSingleMuPrescale3L3PreFiltered" ) Muon2_filterName[_nPair] = 7;
	    	if( _filterName == "hltSingleMuPrescale5L3PreFiltered" ) Muon2_filterName[_nPair] = 8;
	    	if( _filterName == "hltSingleMuPrescale77L3PreFiltered" ) Muon2_filterName[_nPair] = 9;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFiltered9" ) Muon2_filterName[_nPair] = 10;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFiltered11" ) Muon2_filterName[_nPair] = 11;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFiltered13" ) Muon2_filterName[_nPair] = 12;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFiltered15" ) Muon2_filterName[_nPair] = 13;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFiltered" ) Muon2_filterName[_nPair] = 14;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFilteredRelaxedVtx2cm" ) Muon2_filterName[_nPair] = 15;
	    	if( _filterName == "hltSingleMuNoIsoL3PreFilteredRelaxedVtx2mm" ) Muon2_filterName[_nPair] = 16;
	    	if( _filterName == "hltDiMuonIsoL3IsoFiltered" ) Muon2_filterName[_nPair] = 17;
	    	if( _filterName == "hltDiMuonNoIsoL3PreFiltered" ) Muon2_filterName[_nPair] = 18;
	    	if( _filterName == "hltDiMuonNoIsoL3PreFilteredRelaxedVtx2cm" ) Muon2_filterName[_nPair] = 19;
	    	if( _filterName == "hltDiMuonNoIsoL3PreFilteredRelaxedVtx2mm" ) Muon2_filterName[_nPair] = 20;
	    	if( _filterName == "hltJpsiMML3Filtered" ) Muon2_filterName[_nPair] = 21;
	    	if( _filterName == "hltUpsilonMML3Filtered" ) Muon2_filterName[_nPair] = 22;
	    	if( _filterName == "hltZMML3Filtered" ) Muon2_filterName[_nPair] = 23;
	    	if( _filterName == "hltSameSignMuL3PreFiltered" ) Muon2_filterName[_nPair] = 24;
	    	if( _filterName == "hltPsi2SMML3Filtered" ) Muon2_filterName[_nPair] = 25;
	    }
	    */
	    // muon pair index by muon pt
	    MuonPairIndex[_nPair] = _nPair;

	    // highest invariant mass muon pair
	    //bool _isHighestInvMassMuonPair = false;
	    //if( iMuon1->pt() == hMmuon1.pt() && iMuon2->pt() == hMmuon2.pt() ) _isHighestInvMassMuonPair = true;

	    // candidate
	    /* aaa
	    reco::NamedCompositeCandidate aDYcand;
	    aDYcand.setP4(iMuon1->p4() + iMuon2->p4());
	    InvMass[_nPair] = aDYcand.mass();
	    Angle[_nPair] = angleBetween(*iMuon1, *iMuon2);
	    DeltaR[_nPair] = sqrt( pow(iMuon1->eta() - iMuon2->eta(), 2) + pow(iMuon1->phi() - iMuon2->phi(), 2) );
	    Dalpha[_nPair] = (iMuon1->theta() - iMuon2->theta()) + (iMuon1->phi() - iMuon2->phi()) - 2*TMath::Pi();
	    if( iMuon1->charge() != iMuon2->charge() ) isOppSign[_nPair] = 1;
	    else isOppSign[_nPair] = 0;
	    //if( _isHighestInvMassMuonPair ) isHMassPair[_nPair] = 1;
	    aaa */
	    LogDebug("NTUPLE");
	    // muon1 kinematics
	    Muon1_pT[_nPair] = iMuon1->pt();
	    Muon1_Px[_nPair] = iMuon1->px();
	    Muon1_Py[_nPair] = iMuon1->py();
	    Muon1_Pz[_nPair] = iMuon1->pz();
	    Muon1_eta[_nPair] = iMuon1->eta();
	    Muon1_phi[_nPair] = iMuon1->phi();
	    Muon1_sumtrkpt[_nPair] = TrackSumPtrInCone03(*iMuon1, trackHandle);
	    //aaa	    Muon1_trkiso[_nPair] = iMuon1->trackIso();
	    //aaa	    Muon1_hcaliso[_nPair] = iMuon1->hcalIso();
	    //aaa	    Muon1_ecaliso[_nPair] = iMuon1->ecalIso();
	    Muon1_charge[_nPair] = iMuon1->charge();
	    Muon1_nChambers[_nPair] = iMuon1->numberOfChambers(); // # of chambers
	    Muon1_nMatches[_nPair] = iMuon1->numberOfMatches(); // # of chambers with matched segments
	    Muon1_stationMask[_nPair] = iMuon1->stationMask(); // bit map of stations with matched segments
	    // bits 0-1-2-3 = DT stations 1-2-3-4
	    // bits 4-5-6-7 = CSC stations 1-2-3-4
	    int _segments = 0;
	    for( int idet = 1; idet < 4; idet++ ) {
		// DT (1), CSC (2), RPC (3)
		for( int istation = 1; istation < 5; istation++ ) {
		    // station 1, 2, 3, 4
		    _segments += iMuon1->numberOfSegments(istation, idet);
		}
	    }
	    Muon1_nSegments[_nPair] = _segments;
	    Muon1_type[_nPair] = iMuon1->type();
	    if( iMuon1->isGlobalMuon() ) Muon1_muonType[_nPair] = 1; // global muon
	    else if( iMuon1->isStandAloneMuon() ) Muon1_muonType[_nPair] = 2; // STA muon
	    else if( iMuon1->isTrackerMuon() ) Muon1_muonType[_nPair] = 3; // track muon
	    LogDebug("NTUPLE");
	    // global, track and STA
	    reco::TrackRef _muon1Trk;
	    if( Muon1_muonType[_nPair] == 1 ) _muon1Trk = iMuon1->globalTrack();
	    if( Muon1_muonType[_nPair] == 2 ) _muon1Trk = iMuon1->track();
	    if( Muon1_muonType[_nPair] == 3 ) _muon1Trk = iMuon1->standAloneMuon();
	    reco::TrackRef _muon1GlbTrk;
	    reco::TrackRef _muon1StaTrk;
	    reco::TrackRef _muon1TrkTrk;
	    LogDebug("NTUPLE");
	    if(iMuon1->isGlobalMuon()) _muon1GlbTrk = iMuon1->globalTrack();
	    if(iMuon1->isStandAloneMuon()) _muon1StaTrk = iMuon1->outerTrack();
	    if(iMuon1->isTrackerMuon() || iMuon1->isGlobalMuon()) _muon1TrkTrk = iMuon1->innerTrack();
	    LogDebug("NTUPLE");
	    if( _muon1Trk.isNonnull() ) {
	        if( _muon1Trk->normalizedChi2() < 1000 ) Muon1_chi2dof[_nPair] = _muon1Trk->normalizedChi2();
    	        Muon1_nhits[_nPair] = _muon1Trk->numberOfValidHits();
	        Muon1_qoverp[_nPair] = _muon1Trk->qoverp();
	        Muon1_theta[_nPair] = _muon1Trk->theta();
	        Muon1_lambda[_nPair] = _muon1Trk->lambda();
	        Muon1_dxy[_nPair] = _muon1Trk->dxy();
	        Muon1_d0[_nPair] = _muon1Trk->d0();
	        Muon1_dsz[_nPair] = _muon1Trk->dsz();
	        Muon1_dz[_nPair] = _muon1Trk->dz();
	        Muon1_dxyBS[_nPair] = _muon1Trk->dxy(beamSpot.position());
	        Muon1_dszBS[_nPair] = _muon1Trk->dsz(beamSpot.position());
	        Muon1_dzBS[_nPair] = _muon1Trk->dz(beamSpot.position());
	        Muon1_vx[_nPair] = _muon1Trk->vx();
	        Muon1_vy[_nPair] = _muon1Trk->vy();
	        Muon1_vz[_nPair] = _muon1Trk->vz();
	    }

	    if( _muon1GlbTrk.isNonnull() ) {
	        if( _muon1GlbTrk->normalizedChi2() < 1000 ) Muon1_Glb_chi2dof[_nPair] = _muon1GlbTrk->normalizedChi2();
		Muon1_Glb_pT[_nPair] = iMuon1->pt();
		Muon1_Glb_Px[_nPair] = iMuon1->px();
		Muon1_Glb_Py[_nPair] = iMuon1->py();
		Muon1_Glb_Pz[_nPair] = iMuon1->pz();
		Muon1_Glb_eta[_nPair] = iMuon1->eta();
		Muon1_Glb_phi[_nPair] = iMuon1->phi();
		Muon1_Glb_charge[_nPair] = iMuon1->charge();
    	        Muon1_Glb_nhits[_nPair] = _muon1GlbTrk->numberOfValidHits();
	        Muon1_Glb_qoverp[_nPair] = _muon1GlbTrk->qoverp();
	        Muon1_Glb_theta[_nPair] = _muon1GlbTrk->theta();
	        Muon1_Glb_lambda[_nPair] = _muon1GlbTrk->lambda();
	        Muon1_Glb_dxy[_nPair] = _muon1GlbTrk->dxy();
	        Muon1_Glb_d0[_nPair] = _muon1GlbTrk->d0();
	        Muon1_Glb_dsz[_nPair] = _muon1GlbTrk->dsz();
	        Muon1_Glb_dz[_nPair] = _muon1GlbTrk->dz();
	        Muon1_Glb_dxyBS[_nPair] = _muon1GlbTrk->dxy(beamSpot.position());
	        Muon1_Glb_dszBS[_nPair] = _muon1GlbTrk->dsz(beamSpot.position());
	        Muon1_Glb_dzBS[_nPair] = _muon1GlbTrk->dz(beamSpot.position());
	        Muon1_Glb_vx[_nPair] = _muon1GlbTrk->vx();
	        Muon1_Glb_vy[_nPair] = _muon1GlbTrk->vy();
	        Muon1_Glb_vz[_nPair] = _muon1GlbTrk->vz();
		LogDebug("NTUPLE");
		//Now some special stuff
		
		vector<Trajectory> refitted = theRefitter->refit(*_muon1GlbTrk,1) ;

		LogDebug("NTUPLE") << refitted.size();
		std::pair<double,double> thisKink;
		//std::pair<double,double> stachi;
		//std::pair<double,double> tkchi;
		double relative_muon_chi2 = 0.0;
		double relative_tracker_chi2 = 0.0;
		
		if(refitted.size()>0) {
		  thisKink = kink(refitted.front()) ;
		  std::pair<double,double> chi = newChi2(refitted.front());
		  relative_muon_chi2 = chi.second/_muon1StaTrk->ndof();
		  relative_tracker_chi2 = chi.first/_muon1TrkTrk->ndof();
		  
		  LogTrace("NTUPLE") << "thisKink " 
				     << thisKink.first << " " 
				     <<thisKink.second;
		}
		Muon1_Glb_trkKink[_nPair] = thisKink.first;
		Muon1_Glb_glbKink[_nPair] = thisKink.second;
		Muon1_Glb_rChi2Sta[_nPair] = relative_muon_chi2;
		Muon1_Glb_rChi2Trk[_nPair] = relative_tracker_chi2;
		
	    } 

	    if( _muon1StaTrk.isNonnull() ) {
	        if( _muon1StaTrk->normalizedChi2() < 1000 ) Muon1_Sta_chi2dof[_nPair] = _muon1StaTrk->normalizedChi2();
		Muon1_Sta_pT[_nPair] = iMuon1->pt();
		Muon1_Sta_Px[_nPair] = iMuon1->px();
		Muon1_Sta_Py[_nPair] = iMuon1->py();
		Muon1_Sta_Pz[_nPair] = iMuon1->pz();
		Muon1_Sta_eta[_nPair] = iMuon1->eta();
		Muon1_Sta_phi[_nPair] = iMuon1->phi();
		Muon1_Sta_charge[_nPair] = iMuon1->charge();
    	        Muon1_Sta_nhits[_nPair] = _muon1StaTrk->numberOfValidHits();
	        Muon1_Sta_qoverp[_nPair] = _muon1StaTrk->qoverp();
	        Muon1_Sta_theta[_nPair] = _muon1StaTrk->theta();
	        Muon1_Sta_lambda[_nPair] = _muon1StaTrk->lambda();
	        Muon1_Sta_dxy[_nPair] = _muon1StaTrk->dxy();
	        Muon1_Sta_d0[_nPair] = _muon1StaTrk->d0();
	        Muon1_Sta_dsz[_nPair] = _muon1StaTrk->dsz();
	        Muon1_Sta_dz[_nPair] = _muon1StaTrk->dz();
	        Muon1_Sta_dxyBS[_nPair] = _muon1StaTrk->dxy(beamSpot.position());
	        Muon1_Sta_dszBS[_nPair] = _muon1StaTrk->dsz(beamSpot.position());
	        Muon1_Sta_dzBS[_nPair] = _muon1StaTrk->dz(beamSpot.position());
	        Muon1_Sta_vx[_nPair] = _muon1StaTrk->vx();
	        Muon1_Sta_vy[_nPair] = _muon1StaTrk->vy();
	        Muon1_Sta_vz[_nPair] = _muon1StaTrk->vz();
	    }

	    if( _muon1TrkTrk.isNonnull() ) {
	        if( _muon1TrkTrk->normalizedChi2() < 1000 ) Muon1_Trk_chi2dof[_nPair] = _muon1TrkTrk->normalizedChi2();
	    Muon1_Trk_pT[_nPair] = iMuon1->pt();
	    Muon1_Trk_Px[_nPair] = iMuon1->px();
	    Muon1_Trk_Py[_nPair] = iMuon1->py();
	    Muon1_Trk_Pz[_nPair] = iMuon1->pz();
	    Muon1_Trk_eta[_nPair] = iMuon1->eta();
	    Muon1_Trk_phi[_nPair] = iMuon1->phi();
	    Muon1_Trk_charge[_nPair] = iMuon1->charge();
    	        Muon1_Trk_nhits[_nPair] = _muon1TrkTrk->numberOfValidHits();
	        Muon1_Trk_qoverp[_nPair] = _muon1TrkTrk->qoverp();
	        Muon1_Trk_theta[_nPair] = _muon1TrkTrk->theta();
	        Muon1_Trk_lambda[_nPair] = _muon1TrkTrk->lambda();
	        Muon1_Trk_dxy[_nPair] = _muon1TrkTrk->dxy();
	        Muon1_Trk_d0[_nPair] = _muon1TrkTrk->d0();
	        Muon1_Trk_dsz[_nPair] = _muon1TrkTrk->dsz();
	        Muon1_Trk_dz[_nPair] = _muon1TrkTrk->dz();
	        Muon1_Trk_dxyBS[_nPair] = _muon1TrkTrk->dxy(beamSpot.position());
	        Muon1_Trk_dszBS[_nPair] = _muon1TrkTrk->dsz(beamSpot.position());
	        Muon1_Trk_dzBS[_nPair] = _muon1TrkTrk->dz(beamSpot.position());
	        Muon1_Trk_vx[_nPair] = _muon1TrkTrk->vx();
	        Muon1_Trk_vy[_nPair] = _muon1TrkTrk->vy();
	        Muon1_Trk_vz[_nPair] = _muon1TrkTrk->vz();
	    }

	    // MC truthMatch
	    /* aaa
	    if( isMC ) {
            	reco::GenParticleRef genMuon1 = iMuon1->genParticleRef();
		if( genMuon1.isNonnull() ) {
	    	    Muon1_MCtruth_pT[_nPair] = genMuon1->pt();
	    	    Muon1_MCtruth_Px[_nPair] = genMuon1->px();
	    	    Muon1_MCtruth_Py[_nPair] = genMuon1->py();
	    	    Muon1_MCtruth_Pz[_nPair] = genMuon1->pz();
	    	    Muon1_MCtruth_eta[_nPair] = genMuon1->eta();
	    	    Muon1_MCtruth_phi[_nPair] = genMuon1->phi();
	    	    Muon1_MCtruth_charge[_nPair] = genMuon1->charge();
	    	    Muon1_MCtruth_mother[_nPair] = motherId(genMuon1);
		}
	    }
	    aaa */


	    /*
	    // MC truth invariant mass
	    const double par_mass = 0.105658;
	    double genMuon1E = sqrt(Muon1_MCtruth_Px[_nPair]*Muon1_MCtruth_Px[_nPair]
	    				+ Muon1_MCtruth_Py[_nPair]*Muon1_MCtruth_Py[_nPair]
	                                + Muon1_MCtruth_Pz[_nPair]*Muon1_MCtruth_Pz[_nPair]
					+ par_mass*par_mass);
	    TLorentzVector genMuon1vec(Muon1_MCtruth_Px[_nPair], Muon1_MCtruth_Py[_nPair], Muon1_MCtruth_Pz[_nPair], genMuon1E);
	    double genMuon2E = sqrt(Muon2_MCtruth_Px[_nPair]*Muon2_MCtruth_Px[_nPair]
	                                + Muon2_MCtruth_Py[_nPair]*Muon2_MCtruth_Py[_nPair]
					+ Muon2_MCtruth_Pz[_nPair]*Muon2_MCtruth_Pz[_nPair]
					+ par_mass*par_mass);
	    TLorentzVector genMuon2vec(Muon2_MCtruth_Px[_nPair], Muon2_MCtruth_Py[_nPair], Muon2_MCtruth_Pz[_nPair], genMuon2E);
	    TLorentzVector genDYcand = genMuon1vec + genMuon2vec;
	    */
	    _nPair++;
	    //aaa01 }
   }
   // Gen Muon
   if( isMC ) {
       edm::Handle <reco::GenParticleCollection> particles;
       iEvent.getByLabel("genParticles", particles);

       int _GennPair = 0;
       for( size_t m1 = 0; m1 < particles->size(); m1++ ) {
	   const GenParticle &MCmuon1 = (*particles)[m1];
	   if( abs(MCmuon1.pdgId()) != 13 ) continue;
	   //if( fabs(motherId(MCmuon1)) != 23 || fabs(motherId(MCmuon1)) != 22 ) continue;

	   for( size_t m2 = 0; m2 < particles->size(); m2++ ) {
	       const GenParticle &MCmuon2 = (*particles)[m2];

	       if( m1 <= m2 ) continue;
    	       if( abs(MCmuon2.pdgId()) != 13 ) continue;
	       //if( fabs(motherId(MCmuon2)) != 23 || fabs(motherId(MCmuon1)) != 22 ) continue;

		// 4 vectors with muon Gen information
		TVector3 _mu1Gen3V(MCmuon1.px(), MCmuon1.py(), MCmuon1.pz());
		double _mu1mass = MCmuon1.mass();
		double _mu1E = sqrt(_mu1Gen3V.Mag2() + _mu1mass*_mu1mass);
		reco::Particle::LorentzVector _mu1Gen4V(MCmuon1.px(), MCmuon1.py(), MCmuon1.pz(), _mu1E);
	
		TVector3 _mu2Gen3V(MCmuon2.px(), MCmuon2.py(), MCmuon2.pz());
		double _mu2mass = MCmuon2.mass();
		double _mu2E = sqrt(_mu2Gen3V.Mag2() + _mu2mass*_mu2mass);
		reco::Particle::LorentzVector _mu2Gen4V(MCmuon2.px(), MCmuon2.py(), MCmuon2.pz(), _mu2E);

		// candidate
		reco::NamedCompositeCandidate aDYcand;
		aDYcand.setP4(_mu1Gen4V + _mu2Gen4V);
		GENInvMass[_GennPair] = aDYcand.mass();
	
		GENAngle[_GennPair] = angleBetween(MCmuon1, MCmuon2);
		if( MCmuon1.charge() != MCmuon2.charge() ) GENisOppSign[_GennPair] = 1;
		else GENisOppSign[_GennPair] = 0;

		// muon1 kinematics
		GENMuon1_pT[_GennPair] = MCmuon1.pt();
		GENMuon1_Px[_GennPair] = MCmuon1.px();
		GENMuon1_Py[_GennPair] = MCmuon1.py();
		GENMuon1_Pz[_GennPair] = MCmuon1.pz();
		GENMuon1_eta[_GennPair] = MCmuon1.eta();
		GENMuon1_phi[_GennPair] = MCmuon1.phi();
		GENMuon1_charge[_GennPair] = MCmuon1.charge();
		GENMuon1_mother[_GennPair] = motherId(MCmuon1);

		// muon2 kinematics
		GENMuon2_pT[_GennPair] = MCmuon2.pt();
		GENMuon2_Px[_GennPair] = MCmuon2.px();
		GENMuon2_Py[_GennPair] = MCmuon2.py();
		GENMuon2_Pz[_GennPair] = MCmuon2.pz();
		GENMuon2_eta[_GennPair] = MCmuon2.eta();
		GENMuon2_phi[_GennPair] = MCmuon2.phi();
		GENMuon2_charge[_GennPair] = MCmuon2.charge();
		GENMuon2_mother[_GennPair] = motherId(MCmuon2);
	       _GennPair++;
	   }
       }
       GENnPair = _GennPair;
   }
   weight = _intLumi*theFilterEfficiency*theCrossSection/theTotalNevents;
   nPair = _nPair;


   // MET   
   if( metHandle->size() > 1 ) cout << "# of METs = " << metHandle->size() << endl;
   for(edm::View<pat::MET>::const_iterator iMET = metHandle->begin(); iMET != metHandle->end(); ++iMET) {
       MET = iMET->sumEt();
       break;
   }

   // Jets   
   int _njets = 0;
   int _nbjets = 0;
   int _nbjets1 = 0;
   Njets = jetHandle->size();
   for(edm::View<pat::Jet>::const_iterator iJet = jetHandle->begin(); iJet != jetHandle->end(); ++iJet) {
        JETbDiscriminant[_njets] = iJet->bDiscriminator("trackCountingHighEffBJetTags");
	JETflavour[_njets] = iJet->partonFlavour();
	JETcharge[_njets] = iJet->jetCharge();
	JETntracks[_njets] = iJet->associatedTracks().size();
	JETpt[_njets] = iJet->pt();
	JETeta[_njets] = iJet->eta();
	JETphi[_njets] = iJet->phi();

	bool _isbtagged = false;
	if( JETbDiscriminant[_njets] > theBDiscriminant ) {
	    _nbjets++;
	    _isbtagged = true;
	}
	if( _isbtagged ) {
            for(edm::View<reco::Muon>::const_iterator iMuon = muonHandle->begin(); iMuon != muonHandle->end(); ++iMuon) {
		double deltaR = sqrt( pow(iMuon->eta() - iJet->eta(), 2) 
				+ pow(iMuon->phi() - iJet->phi(), 2) );
		if( deltaR < 0.5 ) {
		    _nbjets1++;
		    break;
		}
	    }
	}
	_njets++;
   }
   
   Nbtagged = _nbjets;
   NbtaggedCloseMuon = _nbjets1;
   
   DiMuonTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
AODMuonToNtuple::beginJob(const edm::EventSetup&)
{

  theFile = new TFile(theRootFileName.c_str(),"recreate");
  theFile->cd();

  //
  // DiMuonTree
  //
  DiMuonTree = new TTree("DiMuonTree","DiMuonTree");
  // global event varialbes
  DiMuonTree->Branch("runNum",&runNum,"runNum/I");
  DiMuonTree->Branch("evtNum",&evtNum,"evtNum/I");
  DiMuonTree->Branch("nPair",&nPair,"nPair/I");
  DiMuonTree->Branch("weight",&weight,"weight/D");
  DiMuonTree->Branch("MET", &MET,"MET/D");
  DiMuonTree->Branch("Njets", &Njets,"Njets/I");
  DiMuonTree->Branch("Nmuons", &Nmuons,"Nmuons/I");
  DiMuonTree->Branch("Nbtagged", &Nbtagged,"Nbtagged/I");
  DiMuonTree->Branch("NbtaggedCloseMuon", &NbtaggedCloseMuon,"NbtaggedCloseMuon/I");

  // PV variables
  DiMuonTree->Branch("PVtrackSize", &PVtrackSize,"PVtrackSize/I");
  DiMuonTree->Branch("PVchi2", &PVchi2,"PVchi2/D");
  DiMuonTree->Branch("PVndof", &PVndof,"PVndof/D");
  DiMuonTree->Branch("PVnormalizedChi2", &PVnormalizedChi2,"PVnormalizedChi2/D");
  DiMuonTree->Branch("PVx", &PVx,"PVx/D");
  DiMuonTree->Branch("PVy", &PVy,"PVy/D");
  DiMuonTree->Branch("PVz", &PVz,"PVz/D");

  // trigger list
  DiMuonTree->Branch("HLT_L1Mu", &_HLT_L1Mu,"HLT_L1Mu/I");
  DiMuonTree->Branch("HLT_L1MuOpen", &_HLT_L1MuOpen,"HLT_L1MuOpen/I");
  DiMuonTree->Branch("HLT_L2Mu9", &_HLT_L2Mu9,"HLT_L2Mu9/I");
  DiMuonTree->Branch("HLT_IsoMu9", &_HLT_IsoMu9,"HLT_IsoMu9/I");
  DiMuonTree->Branch("HLT_IsoMu11", &_HLT_IsoMu11,"HLT_IsoMu11/I");
  DiMuonTree->Branch("HLT_IsoMu13", &_HLT_IsoMu13,"HLT_IsoMu13/I");
  DiMuonTree->Branch("HLT_IsoMu15", &_HLT_IsoMu15,"HLT_IsoMu15/I");
  DiMuonTree->Branch("HLT_Mu3", &_HLT_Mu3,"HLT_Mu3/I");
  DiMuonTree->Branch("HLT_Mu5", &_HLT_Mu5,"HLT_Mu5/I");
  DiMuonTree->Branch("HLT_Mu7", &_HLT_Mu7,"HLT_Mu7/I");
  DiMuonTree->Branch("HLT_Mu9", &_HLT_Mu9,"HLT_Mu9/I");
  DiMuonTree->Branch("HLT_Mu11", &_HLT_Mu11,"HLT_Mu11/I");
  DiMuonTree->Branch("HLT_Mu13", &_HLT_Mu13,"HLT_Mu13/I");
  DiMuonTree->Branch("HLT_Mu15", &_HLT_Mu15,"HLT_Mu15/I");
  DiMuonTree->Branch("HLT_Mu15_L1Mu7", &_HLT_Mu15_L1Mu7,"HLT_Mu15_L1Mu7/I");
  DiMuonTree->Branch("HLT_Mu15_Vtx2cm", &_HLT_Mu15_Vtx2cm,"HLT_Mu15_Vtx2cm/I");
  DiMuonTree->Branch("HLT_Mu15_Vtx2mm", &_HLT_Mu15_Vtx2mm,"HLT_Mu15_Vtx2mm/I");
  DiMuonTree->Branch("HLT_DoubleIsoMu3", &_HLT_DoubleIsoMu3,"HLT_DoubleIsoMu3/I");
  DiMuonTree->Branch("HLT_DoubleMu3", &_HLT_DoubleMu3,"HLT_DoubleMu3/I");
  DiMuonTree->Branch("HLT_DoubleMu3_Vtx2cm", &_HLT_DoubleMu3_Vtx2cm,"HLT_DoubleMu3_Vtx2cm/I");
  DiMuonTree->Branch("HLT_DoubleMu3_Vtx2mm", &_HLT_DoubleMu3_Vtx2mm,"HLT_DoubleMu3_Vtx2mm/I");
  DiMuonTree->Branch("HLT_DoubleMu3_JPsi", &_HLT_DoubleMu3_JPsi,"HLT_DoubleMu3_JPsi/I");
  DiMuonTree->Branch("HLT_DoubleMu3_Upsilon", &_HLT_DoubleMu3_Upsilon,"HLT_DoubleMu3_Upsilon/I");
  DiMuonTree->Branch("HLT_DoubleMu7_Z", &_HLT_DoubleMu7_Z,"HLT_DoubleMu7_Z/I");
  DiMuonTree->Branch("HLT_DoubleMu3_SameSign", &_HLT_DoubleMu3_SameSign,"HLT_DoubleMu3_SameSign/I");
  DiMuonTree->Branch("HLT_DoubleMu3_Psi2S", &_HLT_DoubleMu3_Psi2S,"HLT_DoubleMu3_Psi2S/I");

  // Jet
  DiMuonTree->Branch("JETbDiscriminant", &JETbDiscriminant, "JETbDiscriminant[Njets]/D");
  DiMuonTree->Branch("JETflavour", &JETflavour, "JETflavour[Njets]/I");
  DiMuonTree->Branch("JETcharge", &JETcharge, "JETcharge[Njets]/D");
  DiMuonTree->Branch("JETntracks", &JETntracks, "JETntracks[Njets]/I");
  DiMuonTree->Branch("JETpt", &JETpt, "JETpt[Njets]/D");
  DiMuonTree->Branch("JETeta", &JETeta, "JETeta[Njets]/D");
  DiMuonTree->Branch("JETphi", &JETphi, "JETphi[Njets]/D");

  // Global Muon
  // muon pair variables
  DiMuonTree->Branch("InvMass", &InvMass,"InvMass[nPair]/D");
  DiMuonTree->Branch("Angle", &Angle,"Angle[nPair]/D");
  DiMuonTree->Branch("DeltaR", &DeltaR,"DeltaR[nPair]/D");
  DiMuonTree->Branch("Dalpha", &Dalpha,"Dalpha[nPair]/D");
  DiMuonTree->Branch("isOppSign", &isOppSign,"isOppSign[nPair]/I");
  DiMuonTree->Branch("MuonPairIndex", &MuonPairIndex,"MuonPairIndex[nPair]/I");
  DiMuonTree->Branch("isHMassPair", &isHMassPair,"isHMassPair[nPair]/I");
  DiMuonTree->Branch("vtxChi2", &vtxChi2,"vtxChi2[nPair]/D");
  DiMuonTree->Branch("vtxPositionX", &vtxPositionX,"vtxPositionX[nPair]/D");
  DiMuonTree->Branch("vtxPositionY", &vtxPositionY,"vtxPositionY[nPair]/D");
  DiMuonTree->Branch("vtxPositionZ", &vtxPositionZ,"vtxPositionZ[nPair]/D");
  DiMuonTree->Branch("vtxXerror", &vtxXerror,"vtxXerror[nPair]/D");
  DiMuonTree->Branch("vtxYerror", &vtxYerror,"vtxYerror[nPair]/D");
  DiMuonTree->Branch("vtxZerror", &vtxZerror,"vtxZerror[nPair]/D");

  // object variables
  DiMuonTree->Branch("Muon1_type", &Muon1_type,"Muon1_type[nPair]/i");
  DiMuonTree->Branch("Muon1_muonType", &Muon1_muonType,"Muon1_muonType[nPair]/I");
  DiMuonTree->Branch("Muon1_nTrig", &Muon1_nTrig,"Muon1_nTrig[nPair]/I");
  DiMuonTree->Branch("Muon1_triggerObjectType", &Muon1_triggerObjectType,"Muon1_triggerObjectType[nPair]/I");
  DiMuonTree->Branch("Muon1_filterName", &Muon1_filterName,"Muon1_filterName[nPair]/I");
  DiMuonTree->Branch("Muon1_phi", &Muon1_phi,"Muon1_phi[nPair]/D");
  DiMuonTree->Branch("Muon1_eta", &Muon1_eta,"Muon1_eta[nPair]/D");
  DiMuonTree->Branch("Muon1_pT", &Muon1_pT,"Muon1_pT[nPair]/D");
  DiMuonTree->Branch("Muon1_Px", &Muon1_Px,"Muon1_Px[nPair]/D");
  DiMuonTree->Branch("Muon1_Py", &Muon1_Py,"Muon1_Py[nPair]/D");
  DiMuonTree->Branch("Muon1_Pz", &Muon1_Pz,"Muon1_Pz[nPair]/D");
  DiMuonTree->Branch("Muon1_sumtrkpt", &Muon1_sumtrkpt,"Muon1_sumtrkpt[nPair]/D");
  DiMuonTree->Branch("Muon1_trkiso", &Muon1_trkiso,"Muon1_trkiso[nPair]/D");
  DiMuonTree->Branch("Muon1_hcaliso", &Muon1_hcaliso,"Muon1_hcaliso[nPair]/D");
  DiMuonTree->Branch("Muon1_ecaliso", &Muon1_ecaliso,"Muon1_ecaliso[nPair]/D");
  DiMuonTree->Branch("Muon1_charge", &Muon1_charge,"Muon1_charge[nPair]/I");
  DiMuonTree->Branch("Muon1_nChambers", &Muon1_nChambers,"Muon1_nChambers[nPair]/I");
  DiMuonTree->Branch("Muon1_nMatches", &Muon1_nMatches,"Muon1_nMatches[nPair]/I");
  DiMuonTree->Branch("Muon1_stationMask", &Muon1_stationMask,"Muon1_stationMask[nPair]/I");
  DiMuonTree->Branch("Muon1_nSegments", &Muon1_nSegments,"Muon1_nSegments[nPair]/I");
  DiMuonTree->Branch("Muon1_chi2dof", &Muon1_chi2dof,"Muon1_chi2dof[nPair]/D");
  DiMuonTree->Branch("Muon1_nhits", &Muon1_nhits,"Muon1_nhits[nPair]/I");
  DiMuonTree->Branch("Muon1_qoverp", &Muon1_qoverp,"Muon1_qoverp[nPair]/D");
  DiMuonTree->Branch("Muon1_theta", &Muon1_theta,"Muon1_theta[nPair]/D");
  DiMuonTree->Branch("Muon1_lambda", &Muon1_lambda,"Muon1_lambda[nPair]/D");
  DiMuonTree->Branch("Muon1_dxy", &Muon1_dxy,"Muon1_dxy[nPair]/D");
  DiMuonTree->Branch("Muon1_d0", &Muon1_d0,"Muon1_d0[nPair]/D");
  DiMuonTree->Branch("Muon1_dsz", &Muon1_dsz,"Muon1_dsz[nPair]/D");
  DiMuonTree->Branch("Muon1_dz", &Muon1_dz,"Muon1_dz[nPair]/D");
  DiMuonTree->Branch("Muon1_dxyBS", &Muon1_dxyBS,"Muon1_dxyBS[nPair]/D");
  DiMuonTree->Branch("Muon1_dszBS", &Muon1_dszBS,"Muon1_dszBS[nPair]/D");
  DiMuonTree->Branch("Muon1_dzBS", &Muon1_dzBS,"Muon1_dzBS[nPair]/D");
  DiMuonTree->Branch("Muon1_vx", &Muon1_vx,"Muon1_vx[nPair]/D");
  DiMuonTree->Branch("Muon1_vy", &Muon1_vy,"Muon1_vy[nPair]/D");
  DiMuonTree->Branch("Muon1_vz", &Muon1_vz,"Muon1_vz[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_chi2dof", &Muon1_Glb_chi2dof,"Muon1_Glb_chi2dof[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_phi", &Muon1_Glb_phi,"Muon1_Glb_phi[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_eta", &Muon1_Glb_eta,"Muon1_Glb_eta[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_pT", &Muon1_Glb_pT,"Muon1_Glb_pT[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_Px", &Muon1_Glb_Px,"Muon1_Glb_Px[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_Py", &Muon1_Glb_Py,"Muon1_Glb_Py[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_Pz", &Muon1_Glb_Pz,"Muon1_Glb_Pz[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_charge", &Muon1_Glb_charge,"Muon1_Glb_charge[nPair]/I");
  DiMuonTree->Branch("Muon1_Glb_nhits", &Muon1_Glb_nhits,"Muon1_Glb_nhits[nPair]/I");
  DiMuonTree->Branch("Muon1_Glb_qoverp", &Muon1_Glb_qoverp,"Muon1_Glb_qoverp[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_theta", &Muon1_Glb_theta,"Muon1_Glb_theta[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_lambda", &Muon1_Glb_lambda,"Muon1_Glb_lambda[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_dxy", &Muon1_Glb_dxy,"Muon1_Glb_dxy[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_d0", &Muon1_Glb_d0,"Muon1_Glb_d0[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_dsz", &Muon1_Glb_dsz,"Muon1_Glb_dsz[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_dz", &Muon1_Glb_dz,"Muon1_Glb_dz[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_dxyBS", &Muon1_Glb_dxyBS,"Muon1_Glb_dxyBS[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_dszBS", &Muon1_Glb_dszBS,"Muon1_Glb_dszBS[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_dzBS", &Muon1_Glb_dzBS,"Muon1_Glb_dzBS[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_vx", &Muon1_Glb_vx,"Muon1_Glb_vx[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_vy", &Muon1_Glb_vy,"Muon1_Glb_vy[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_vz", &Muon1_Glb_vz,"Muon1_Glb_vz[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_trkKink", &Muon1_Glb_trkKink,"Muon1_Glb_trkKink[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_glbKink", &Muon1_Glb_glbKink,"Muon1_Glb_glbKink[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_rChi2Sta", &Muon1_Glb_rChi2Sta,"Muon1_Glb_rChi2Sta[nPair]/D");
  DiMuonTree->Branch("Muon1_Glb_rChi2Trk", &Muon1_Glb_rChi2Trk,"Muon1_Glb_rChi2Trk[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_chi2dof", &Muon1_Sta_chi2dof,"Muon1_Sta_chi2dof[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_phi", &Muon1_Sta_phi,"Muon1_Sta_phi[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_eta", &Muon1_Sta_eta,"Muon1_Sta_eta[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_pT", &Muon1_Sta_pT,"Muon1_Sta_pT[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_Px", &Muon1_Sta_Px,"Muon1_Sta_Px[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_Py", &Muon1_Sta_Py,"Muon1_Sta_Py[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_Pz", &Muon1_Sta_Pz,"Muon1_Sta_Pz[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_charge", &Muon1_Sta_charge,"Muon1_Sta_charge[nPair]/I");
  DiMuonTree->Branch("Muon1_Sta_nhits", &Muon1_Sta_nhits,"Muon1_Sta_nhits[nPair]/I");
  DiMuonTree->Branch("Muon1_Sta_qoverp", &Muon1_Sta_qoverp,"Muon1_Sta_qoverp[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_theta", &Muon1_Sta_theta,"Muon1_Sta_theta[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_lambda", &Muon1_Sta_lambda,"Muon1_Sta_lambda[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_dxy", &Muon1_Sta_dxy,"Muon1_Sta_dxy[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_d0", &Muon1_Sta_d0,"Muon1_Sta_d0[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_dsz", &Muon1_Sta_dsz,"Muon1_Sta_dsz[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_dz", &Muon1_Sta_dz,"Muon1_Sta_dz[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_dxyBS", &Muon1_Sta_dxyBS,"Muon1_Sta_dxyBS[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_dszBS", &Muon1_Sta_dszBS,"Muon1_Sta_dszBS[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_dzBS", &Muon1_Sta_dzBS,"Muon1_Sta_dzBS[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_vx", &Muon1_Sta_vx,"Muon1_Sta_vx[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_vy", &Muon1_Sta_vy,"Muon1_Sta_vy[nPair]/D");
  DiMuonTree->Branch("Muon1_Sta_vz", &Muon1_Sta_vz,"Muon1_Sta_vz[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_chi2dof", &Muon1_Trk_chi2dof,"Muon1_Trk_chi2dof[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_phi", &Muon1_Trk_phi,"Muon1_Trk_phi[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_eta", &Muon1_Trk_eta,"Muon1_Trk_eta[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_pT", &Muon1_Trk_pT,"Muon1_Trk_pT[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_Px", &Muon1_Trk_Px,"Muon1_Trk_Px[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_Py", &Muon1_Trk_Py,"Muon1_Trk_Py[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_Pz", &Muon1_Trk_Pz,"Muon1_Trk_Pz[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_charge", &Muon1_Trk_charge,"Muon1_Trk_charge[nPair]/I");
  DiMuonTree->Branch("Muon1_Trk_nhits", &Muon1_Trk_nhits,"Muon1_Trk_nhits[nPair]/I");
  DiMuonTree->Branch("Muon1_Trk_qoverp", &Muon1_Trk_qoverp,"Muon1_Trk_qoverp[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_theta", &Muon1_Trk_theta,"Muon1_Trk_theta[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_lambda", &Muon1_Trk_lambda,"Muon1_Trk_lambda[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_dxy", &Muon1_Trk_dxy,"Muon1_Trk_dxy[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_d0", &Muon1_Trk_d0,"Muon1_Trk_d0[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_dsz", &Muon1_Trk_dsz,"Muon1_Trk_dsz[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_dz", &Muon1_Trk_dz,"Muon1_Trk_dz[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_dxyBS", &Muon1_Trk_dxyBS,"Muon1_Trk_dxyBS[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_dszBS", &Muon1_Trk_dszBS,"Muon1_Trk_dszBS[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_dzBS", &Muon1_Trk_dzBS,"Muon1_Trk_dzBS[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_vx", &Muon1_Trk_vx,"Muon1_Trk_vx[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_vy", &Muon1_Trk_vy,"Muon1_Trk_vy[nPair]/D");
  DiMuonTree->Branch("Muon1_Trk_vz", &Muon1_Trk_vz,"Muon1_Trk_vz[nPair]/D");
  /*
  DiMuonTree->Branch("Muon1_MCtruth_pT", &Muon1_MCtruth_pT,"Muon1_MCtruth_pT[nPair]/D");
  DiMuonTree->Branch("Muon1_MCtruth_Px", &Muon1_MCtruth_Px,"Muon1_MCtruth_Px[nPair]/D");
  DiMuonTree->Branch("Muon1_MCtruth_Py", &Muon1_MCtruth_Py,"Muon1_MCtruth_Py[nPair]/D");
  DiMuonTree->Branch("Muon1_MCtruth_Pz", &Muon1_MCtruth_Pz,"Muon1_MCtruth_Pz[nPair]/D");
  DiMuonTree->Branch("Muon1_MCtruth_eta", &Muon1_MCtruth_eta,"Muon1_MCtruth_eta[nPair]/D");
  DiMuonTree->Branch("Muon1_MCtruth_phi", &Muon1_MCtruth_phi,"Muon1_MCtruth_phi[nPair]/D");
  DiMuonTree->Branch("Muon1_MCtruth_charge", &Muon1_MCtruth_charge,"Muon1_MCtruth_charge[nPair]/I");
  DiMuonTree->Branch("Muon1_MCtruth_mother", &Muon1_MCtruth_mother,"Muon1_MCtruth_mother[nPair]/I");
  */

  // Gen Muon
  // global event varialbes (for GEN)
  DiMuonTree->Branch("GENnPair",&GENnPair,"GENnPair/I");
  // muon pair variables
  DiMuonTree->Branch("GENInvMass", &GENInvMass,"GENInvMass[GENnPair]/D");
  DiMuonTree->Branch("GENAngle", &GENAngle,"GENAngle[GENnPair]/D");
  DiMuonTree->Branch("GENisOppSign", &GENisOppSign,"GENisOppSign[GENnPair]/I");
  // object variables
  DiMuonTree->Branch("GENMuon1_phi", &GENMuon1_phi,"GENMuon1_phi[GENnPair]/D");
  DiMuonTree->Branch("GENMuon1_eta", &GENMuon1_eta,"GENMuon1_eta[GENnPair]/D");
  DiMuonTree->Branch("GENMuon1_pT", &GENMuon1_pT,"GENMuon1_pT[GENnPair]/D");
  DiMuonTree->Branch("GENMuon1_Px", &GENMuon1_Px,"GENMuon1_Px[GENnPair]/D");
  DiMuonTree->Branch("GENMuon1_Py", &GENMuon1_Py,"GENMuon1_Py[GENnPair]/D");
  DiMuonTree->Branch("GENMuon1_Pz", &GENMuon1_Pz,"GENMuon1_Pz[GENnPair]/D");
  DiMuonTree->Branch("GENMuon1_mother", &GENMuon1_mother,"GENMuon1_mother[GENnPair]/D");
  DiMuonTree->Branch("GENMuon1_charge", &GENMuon1_charge,"GENMuon1_charge[GENnPair]/I");
  DiMuonTree->Branch("GENMuon2_phi", &GENMuon2_phi,"GENMuon2_phi[GENnPair]/D");
  DiMuonTree->Branch("GENMuon2_eta", &GENMuon2_eta,"GENMuon2_eta[GENnPair]/D");
  DiMuonTree->Branch("GENMuon2_pT", &GENMuon2_pT,"GENMuon2_pT[GENnPair]/D");
  DiMuonTree->Branch("GENMuon2_Px", &GENMuon2_Px,"GENMuon2_Px[GENnPair]/D");
  DiMuonTree->Branch("GENMuon2_Py", &GENMuon2_Py,"GENMuon2_Py[GENnPair]/D");
  DiMuonTree->Branch("GENMuon2_Pz", &GENMuon2_Pz,"GENMuon2_Pz[GENnPair]/D");
  DiMuonTree->Branch("GENMuon2_mother", &GENMuon2_mother,"GENMuon2_mother[GENnPair]/D");
  DiMuonTree->Branch("GENMuon2_charge", &GENMuon2_charge,"GENMuon2_charge[GENnPair]/I");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AODMuonToNtuple::endJob() {
  std::cout <<"++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout <<"analyzed " << nEvt << " events: " << std::endl;
  std::cout << "writing information into file: " << theFile->GetName() << std::endl;
  std::cout <<"++++++++++++++++++++++++++++++++++++++" << std::endl;

  theFile->Write();
  theFile->Close();

}

float AODMuonToNtuple::deltaEtaPhi(const reco::Muon& lhs, const reco::Muon& rhs) const {
    float phi1 = lhs.track()->phi();
    float eta1 = lhs.track()->eta();

    float phi2 = rhs.track()->phi();
    float eta2 = rhs.track()->eta();
    return deltaR(eta1, phi1, eta2, phi2);

}

float AODMuonToNtuple::angleBetween(const reco::Muon& lhs, const reco::Muon& rhs) const {
   GlobalVector mom1(lhs.track()->px(), lhs.track()->py(), lhs.track()->pz());
   GlobalVector mom2(rhs.track()->px(), rhs.track()->py(), rhs.track()->pz());

   GlobalVector dmom = mom1 - mom2;
   return acos( ( mom1.mag() * mom1.mag() + mom2.mag() * mom2.mag() - dmom.mag() * dmom.mag() ) / (2*mom1.mag()*mom2.mag() ));

}

float AODMuonToNtuple::angleBetween(const reco::Track& lhs, const reco::Track& rhs) const {
   GlobalVector mom1(lhs.px(), lhs.py(), lhs.pz());
   GlobalVector mom2(rhs.px(), rhs.py(), rhs.pz());

   GlobalVector dmom = mom1 - mom2;
   return acos( ( mom1.mag() * mom1.mag() + mom2.mag() * mom2.mag() - dmom.mag() * dmom.mag() ) / (2*mom1.mag()*mom2.mag() ));

}

float AODMuonToNtuple::angleBetween(const reco::GenParticle& lhs, const reco::GenParticle& rhs) const {
   GlobalVector mom1(lhs.px(), lhs.py(), lhs.pz());
   GlobalVector mom2(rhs.px(), rhs.py(), rhs.pz());

   GlobalVector dmom = mom1 - mom2;
   return acos( ( mom1.mag() * mom1.mag() + mom2.mag() * mom2.mag() - dmom.mag() * dmom.mag() ) / (2*mom1.mag()*mom2.mag() ));

}

int AODMuonToNtuple::motherId(const reco::GenParticle& par) const
{
    int _motherid = -9999;
    const Candidate* motherCand = par.mother(0);

    if( !motherCand ) return -9999;
    while( motherCand ) {
	if(  fabs(motherCand->pdgId()) == 2212 ) break;
	if(  fabs(motherCand->pdgId()) < 7 ) break;
	_motherid = motherCand->pdgId();
	motherCand = motherCand->mother(0);
	if( fabs(_motherid) == 23 || fabs(_motherid) == 22 ) return _motherid;
    }
    return _motherid;
}

int AODMuonToNtuple::motherId(const reco::GenParticleRef par) const
{
    int _motherid = -9999;
    const Candidate* motherCand = par->mother(0);

    if( !motherCand ) return -9999;
    while( motherCand ) {
	if(  fabs(motherCand->pdgId()) == 2212 ) break;
	if(  fabs(motherCand->pdgId()) < 7 ) break;
	_motherid = motherCand->pdgId();
	motherCand = motherCand->mother(0);
	if( fabs(_motherid) == 23 || fabs(_motherid) == 22 ) return _motherid;
    }
    return _motherid;
}

double AODMuonToNtuple::TrackSumPtrInCone03( const reco::Muon& muon, edm::Handle<edm::View<reco::Track> > tracks )
{
   double sum_pt = 0;
   for(edm::View<reco::Track>::const_iterator itrk = tracks->begin(); itrk != tracks->end(); ++itrk) {
       double dR = deltaR(muon.eta(), itrk->eta(), muon.phi(), itrk->phi());
       if( dR < 0.3 ) {
	   sum_pt += itrk->pt();
       }
   }
   return sum_pt;
}

std::pair<double,double> AODMuonToNtuple::newChi2(Trajectory& muon) const {
  double muChi2 = 0.0;
  double tkChi2 = 0.0;
  
  typedef TransientTrackingRecHit::ConstRecHitPointer 	ConstRecHitPointer;
  typedef ConstRecHitPointer RecHit;
  typedef vector<TrajectoryMeasurement>::const_iterator TMI;

  vector<TrajectoryMeasurement> meas = muon.measurements();

  for ( TMI m = meas.begin(); m != meas.end(); m++ ) {
    TransientTrackingRecHit::ConstRecHitPointer hit = m->recHit();
    const TrajectoryStateOnSurface& uptsos = (*m).updatedState();
    TransientTrackingRecHit::RecHitPointer preciseHit = hit->clone(uptsos);
    double estimate = 0.0;
    if (preciseHit->isValid() && uptsos.isValid()) {
      estimate = theEstimator->estimate(uptsos, *preciseHit ).second;
    }
    
    LogTrace("NTUPLE") << "estimate " << estimate << " TM.est " << m->estimate();
    double tkDiff = 0.0;
    double staDiff = 0.0;
    if ( hit->isValid() &&  (hit->geographicalId().det()) == DetId::Tracker ) {
      tkChi2 += estimate;
      tkDiff = estimate - m->estimate();
    }
    if ( hit->isValid() &&  (hit->geographicalId().det()) == DetId::Muon ) {
      muChi2 += estimate;
      staDiff = estimate - m->estimate();
    }
  }
  
  return std::pair<double,double>(tkChi2,muChi2);
       
}

std::pair<double,double> AODMuonToNtuple::kink(Trajectory& muon) const {
    
  double result = 0.0;
  double resultGlb = 0.0;
  
  typedef TransientTrackingRecHit::ConstRecHitPointer 	ConstRecHitPointer;
  typedef ConstRecHitPointer RecHit;
  typedef vector<TrajectoryMeasurement>::const_iterator TMI;

  vector<TrajectoryMeasurement> meas = muon.measurements();
  //LogDebug("NTUPLE");
  
  for ( TMI m = meas.begin(); m != meas.end(); m++ ) {
    //LogDebug("NTUPLE");

    TransientTrackingRecHit::ConstRecHitPointer hit = m->recHit();

    //double estimate = 0.0;

    RecHit rhit = (*m).recHit();
    bool ok = false;
    if ( rhit->isValid() ) {
      if(DetId::Tracker == rhit->geographicalId().det()) ok = true;
    }

    //if ( !ok ) continue;
    
    const TrajectoryStateOnSurface& tsos = (*m).predictedState();


    if ( tsos.isValid() ) {

      double phi1 = tsos.globalPosition().phi();
      if ( phi1 < 0 ) phi1 = 2*M_PI + phi1;

      double phi2 = rhit->globalPosition().phi();
      if ( phi2 < 0 ) phi2 = 2*M_PI + phi2;

      double diff = fabs(phi1 - phi2);
      if ( diff > M_PI ) diff = 2*M_PI - diff;

      GlobalPoint hitPos = rhit->globalPosition();

      GlobalError hitErr = rhit->globalPositionError();

      double error = hitErr.phierr(hitPos);  // error squared

      double s = ( error > 0.0 ) ? (diff*diff)/error : (diff*diff);

      if(ok) result += s;
      resultGlb += s;
    }   
  }
  
  return pair<double,double>(result,resultGlb);

}

//define this as a plug-in
//DEFINE_FWK_MODULE(AODMuonToNtuple);