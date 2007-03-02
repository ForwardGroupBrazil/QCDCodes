// -*- C++ -*-
//
// Package:    TrackingRegionAnalyzer
// Class:      TrackingRegionAnalyzer
// 
/**\class TrackingRegionAnalyzer TrackingRegionAnalyzer.cc UserCode/TrackingRegionAnalyzer/src/TrackingRegionAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam Everett
//         Created:  Tue Feb 20 17:29:48 CET 2007
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

//Data Formats
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Common/interface/Ref.h"

//used object classes
//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
//#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "RecoTracker/TkTrackingRegions/interface/RectangularEtaPhiTrackingRegion.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "TrackingTools/GeomPropagators/interface/StateOnTrackerBound.h"

#include <TH1.h>
#include <TH2.h>
#include <vector>
#include <algorithm>

//
// class decleration
//

class TrackingRegionAnalyzer : public edm::EDAnalyzer {
public:
  explicit TrackingRegionAnalyzer(const edm::ParameterSet&);
  ~TrackingRegionAnalyzer();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  RectangularEtaPhiTrackingRegion defineRegionOfInterest1(const reco::TrackRef&) const ;
  RectangularEtaPhiTrackingRegion defineRegionOfInterest2(const reco::TrackRef&) const ;
  RectangularEtaPhiTrackingRegion defineRegionOfInterest3(const reco::TrackRef&) const ;
  reco::TrackCollection chooseRegionalTrackerTracks(const RectangularEtaPhiTrackingRegion&, const reco::TrackCollection&) const ;

  // ----------member data ---------------------------
  MuonServiceProxy* theService;
  std::string stateOnTrackerOutProp;

  edm::InputTag theSTALabel;
  edm::InputTag theGLBLabel;
  edm::InputTag theTKLabel;
  edm::InputTag theSIMLabel;

  TH2F *h_region1, *h_region2, *h_region3;
  TH2F *h_region1a, *h_region2a, *h_region3a, *h_region3b;
  TH1 *h_tracks0, *h_tracks1, *h_tracks2, *h_tracks3, *h_n_sta0, *h_n_sta1, *h_n_tk1;
  TH1 *h_optimal, *h_oldRegion;
  TH2 *h_etaFactor, *h_phiFactor;

  TH1 *h_deta, *h_dphi;
  TH2 *h_deta1, *h_deta2, *h_dphi1, *h_dphi2, *h_deta3, *h_dphi3;
  TH2 *h_deta4, *h_deta5, *h_dphi4, *h_dphi5, *h_deta6, *h_dphi6;

};

using namespace edm;
using namespace reco;
using namespace std;

//
// constants, enums and typedefs
//
typedef PixelRecoRange< float > Range;
typedef TkTrackingRegionsMargin< float > Margin;

//
// static data member definitions
//

//
// constructors and destructor
//
TrackingRegionAnalyzer::TrackingRegionAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  theSTALabel = iConfig.getParameter<edm::InputTag>("STACollectionLabel");
  theGLBLabel = iConfig.getParameter<edm::InputTag>("GLBCollectionLabel");
  theTKLabel = iConfig.getParameter<edm::InputTag>("TKCollectionLabel");
  theSIMLabel = iConfig.getParameter<edm::InputTag>("SIMCollectionLabel");

  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  LogDebug("Analyzer");  

  stateOnTrackerOutProp = iConfig.getParameter<std::string>("StateOnTrackerBoundOutPropagator");
}


TrackingRegionAnalyzer::~TrackingRegionAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  if(theService) delete theService;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TrackingRegionAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  const string category = "Analyzer|TrackingRegionAnalyzer";
  using namespace edm;

  LogDebug(category);
  
  theService->update(iSetup);

 Handle<edm::SimTrackContainer> SIMTrackCollection;
 iEvent.getByLabel(theSIMLabel,SIMTrackCollection);
 const SimTrackContainer simTC = *(SIMTrackCollection.product());
 
 Handle<reco::TrackCollection> TKTrackCollection;
 iEvent.getByLabel( theTKLabel, TKTrackCollection);
 const reco::TrackCollection tkTC = *(TKTrackCollection.product());
 
 Handle<reco::TrackCollection> STATrackCollection;
 iEvent.getByLabel( theSTALabel, STATrackCollection);
 const reco::TrackCollection staTC = *(STATrackCollection.product());

 Handle<reco::TrackCollection> GLBTrackCollection;
 iEvent.getByLabel( theGLBLabel, GLBTrackCollection);
 const reco::TrackCollection glbTC = *(GLBTrackCollection.product());
 
 Handle<reco::MuonCollection> MuCollection;
 iEvent.getByLabel(theGLBLabel,MuCollection);
 const reco::MuonCollection muonC = *(MuCollection.product());

 h_tracks0->Fill(tkTC.size());
 h_n_sta0->Fill(staTC.size());

 for(MuonCollection::const_iterator mu = muonC.begin(); mu != muonC.end(); ++mu) {
   TrackRef sta = mu->standAloneMuon();
   TrackRef trk = mu->track();

   float deta = (fabs(trk->eta() - sta->eta()));
   float dphi = (fabs(Geom::Phi<float>(trk->phi())-Geom::Phi<float>(sta->phi())));

   if(trk->eta() > 2.5 || sta->eta() > 2.5) LogDebug(category) << "*********\n**********\n**********\nEta out of range trk: " << trk->eta() << " sta: " << sta->eta() << "\n**********\n**********\n**********";

   h_deta->Fill(deta);
   h_dphi->Fill(dphi);
   h_optimal->Fill(deta,dphi);
   h_deta1->Fill(fabs(sta->eta()),deta);
   h_dphi1->Fill(fabs(sta->eta()),dphi);

   h_deta2->Fill(sta->pt(),deta);
   h_dphi2->Fill(sta->pt(),dphi);

   RectangularEtaPhiTrackingRegion regionX = defineRegionOfInterest3(sta);

   Range  etaRange  = regionX.etaRange();
   Margin phiMargin = regionX.phiMargin();
   
   float regiondeta =  fabs(etaRange.max()-etaRange.mean());
   float regiondphi =  phiMargin.right();

   h_deta3->Fill(deta,regiondeta);
   h_dphi3->Fill(dphi,regiondphi);

   h_deta4->Fill(fabs(sta->eta()),regiondeta);
   h_dphi4->Fill(fabs(sta->eta()),regiondphi);

   h_deta5->Fill(sta->pt(),regiondeta);
   h_dphi5->Fill(sta->pt(),regiondphi);
   
   TrackCollection regionalTkTracks =  chooseRegionalTrackerTracks(regionX,tkTC);

   RectangularEtaPhiTrackingRegion region1 = defineRegionOfInterest1(sta);
   TrackCollection regionalTkTracks1 =  chooseRegionalTrackerTracks(region1,tkTC);
   RectangularEtaPhiTrackingRegion region2 = defineRegionOfInterest2(sta);
   TrackCollection regionalTkTracks2 =  chooseRegionalTrackerTracks(region2,tkTC);
   LogInfo(category) << "TkTracks in region: " <<regionalTkTracks.size();
   h_tracks1->Fill(regionalTkTracks1.size());    
   h_tracks2->Fill(regionalTkTracks2.size());    
   h_tracks3->Fill(regionalTkTracks.size());    

   if(regionalTkTracks.size() < 1 ) {
     LogDebug(category) << "&&&&&&& " << regionalTkTracks1.size() 
			<< " " << regionalTkTracks2.size()
			<< " " << regionalTkTracks.size();
     LogDebug(category) << "&&&&&&& " << deta << " " << dphi;
     LogDebug(category) << "&&&&&&& Trk" << trk->eta() << " " 
			<< trk->phi() << " " << trk->pt();
     LogDebug(category) << "&&&&&&& Sta" << sta->eta() << " " 
			<< sta->phi() << " " << sta->pt();
   }
   
 }
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackingRegionAnalyzer::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;

  h_deta =  fs->make<TH1F>("deta","#Delta #eta",300,0,3);
  h_dphi =  fs->make<TH1F>("dphi","#Delta #phi",300,0,3);

  h_deta1 =  fs->make<TH2F>("deta1","#Delta #eta vs #eta",30,0,3,300,0,3);
  h_dphi1 =  fs->make<TH2F>("dphi1","#Delta #phi vs #eta",30,0,3,300,0,3);

  h_deta2 =  fs->make<TH2F>("deta2","#Delta #eta vs p_{t}",20,0,200,300,0,3);
  h_dphi2 =  fs->make<TH2F>("dphi2","#Delta #phi vs p_{t}",20,0,200,300,0,3);

  h_deta3 =  fs->make<TH2F>("deta3","#sigma_{#eta} vs #Delta #eta",300,0,3,300,0,3);
  h_dphi3 =  fs->make<TH2F>("dphi3","#sigma_{#phi} vs #Delta #phi",300,0,3,300,0,3);

  h_deta4 =  fs->make<TH2F>("deta4","#sigma_{#eta} vs #eta",30,0,3,300,0,3);
  h_dphi4 =  fs->make<TH2F>("dphi4","#sigma_{#phi} vs #eta",30,0,3,300,0,3);

  h_deta5 =  fs->make<TH2F>("deta5","#sigma_{#eta} vs p_{t}",20,0,200,300,0,3);
  h_dphi5 =  fs->make<TH2F>("dphi5","#sigma_{#phi} vs p_{t}",20,0,200,300,0,3);

  h_region1 =  fs->make<TH2F>("regionSize1","Size of Region",30,0,3,30,0,3);
  h_region2 =  fs->make<TH2F>("regionSize2","Size of Region",30,0,3,30,0,3);
  h_region3 =  fs->make<TH2F>("regionSize3","Size of Region",30,0,3,30,0,3);

  h_region1a =  fs->make<TH2F>("regionSize1a","Size of Region",300,0,3,300,0,3);
  h_region2a =  fs->make<TH2F>("regionSize2a","Size of Region",300,0,3,300,0,3);
  h_region3a =  fs->make<TH2F>("regionSize3a","Size of Region",300,0,3,300,0,3);
  h_region3b =  fs->make<TH2F>("regionSize3b","Size of Region",300,0,3,300,0,3);

  h_optimal =  fs->make<TH2F>("optimal","Optimal Size of Region",30,0,3,30,0,3);
  h_oldRegion =  fs->make<TH2F>("oldRegion","Old Size of Region",30,0,3,30,0,3);
  
  h_etaFactor =  fs->make<TH2F>("etaFactor","Size of Region",30,0,3,30,0,30);
  h_phiFactor =  fs->make<TH2F>("phiFactor","Size of Region",30,0,3,30,0,30);

  h_tracks0 = fs->make<TH1F>("tracks0","Tk Tracks in Event",21,-0.5,20.5);
  h_tracks1 = fs->make<TH1F>("tracks1","Tk Tracks in Event",21,-0.5,20.5);
  h_tracks2 = fs->make<TH1F>("tracks2","Tk Tracks in Event",21,-0.5,20.5);
  h_tracks3 = fs->make<TH1F>("tracks3","Tk Tracks in Event",21,-0.5,20.5);

  h_n_sta0 = fs->make<TH1F>("n_sta0","STA Tracks in Event",21,-0.5,20.5);
  h_n_sta1 = fs->make<TH1F>("n_sta1","STA Tracks in Event",21,-0.5,20.5);
  h_n_tk1 = fs->make<TH1F>("n_tk1","Tk Tracks in Event",21,-0.5,20.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackingRegionAnalyzer::endJob() {
}

//
// define a region of interest within the tracker
//
RectangularEtaPhiTrackingRegion TrackingRegionAnalyzer::defineRegionOfInterest1(const reco::TrackRef& staTrack) const {
  //region 1 is FTS momentum - TSOS position -- with mods

  TrajectoryStateTransform tsTransform;
  FreeTrajectoryState muFTS = tsTransform.initialFreeState(*staTrack,&*theService->magneticField());
  
  //Get Track direction at vertex
  GlobalVector dirVector(staTrack->px(),staTrack->py(),staTrack->pz());
  LogDebug("Analyzer") << dirVector;
  //--
  dirVector = muFTS.momentum();
  LogDebug("Analyzer") << dirVector;
  
  //Get track momentum
  const math::XYZVector& mo = staTrack->innerMomentum();
  GlobalVector mom(mo.x(),mo.y(),mo.z());
  if ( staTrack->p() > 1.0 ) {
    mom = dirVector; 
  }
  //Get innerMu position
  const math::XYZPoint& po = staTrack->innerPosition();
  GlobalPoint pos(po.x(),po.y(),po.z());
  
  //Get dEta and dPhi: (direction at vertex) - (innerMuTsos position)
  float eta1 = dirVector.eta();
  float eta2 = pos.eta();
  float deta(fabs(eta1- eta2));
  float dphi(fabs(Geom::Phi<float>(dirVector.phi())-Geom::Phi<float>(pos.phi())));
  
  //deta = 1 * deta;
  //dphi = 1 * dphi;
  
  //deta = 1 * max(double(deta),0.05);
  //dphi = 1 * max(double(dphi),0.07);
  
  GlobalPoint vertexPos = (muFTS.position());
  GlobalError vertexErr = (muFTS.cartesianError().position());
  
  double minPt    = max(1.5,mom.perp()*0.6);
  double deltaZ   = min(15.9,3*sqrt(vertexErr.czz()));
  double deltaEta = 0.05;
  double deltaPhi = 0.07;
  
  if ( deta > 0.05 ) { // 0.06
    deltaEta += deta/2;
  }
  if ( dphi > 0.07 ) {
    deltaPhi += 0.15;
    if ( fabs(eta2) < 1.0 && mom.perp() < 6. ) deltaPhi = dphi;
  }
  if ( fabs(eta1) < 1.25 && fabs(eta1) > 0.8 ) deltaEta = max(0.07,deltaEta);
  if ( fabs(eta1) < 1.3  && fabs(eta1) > 1.0 ) deltaPhi = max(0.3,deltaPhi);

  deltaEta = 1 * max(double(2.5 * deta),deltaEta);
  deltaPhi = 1 * max(double(3.5 * dphi),deltaPhi);
  //deltaEta = deta;
  //deltaPhi = dphi;

  RectangularEtaPhiTrackingRegion rectRegion(dirVector, vertexPos,
                                             minPt, 0.2,
                                             deltaZ, deltaEta, deltaPhi);
  h_region1->Fill(deltaEta,deltaPhi);
  h_region1a->Fill(deta,dphi);
  return rectRegion;

}

//
// define a region of interest within the tracker
//
RectangularEtaPhiTrackingRegion TrackingRegionAnalyzer::defineRegionOfInterest2(const reco::TrackRef& staTrack) const {

  //region 2 is FTS mometum and TSOS position deltaEta deltaPhi

  // Get Mu FTS updated at vertex
  TrajectoryStateTransform tsTransform;
  FreeTrajectoryState muFTS = tsTransform.initialFreeState(*staTrack,&*theService->magneticField());
  
  //Get Track direction at vertex
  GlobalVector dirVector(muFTS.momentum());
  
  //Get track momentum
  const math::XYZVector& mo = staTrack->innerMomentum();
  GlobalVector mom(mo.x(),mo.y(),mo.z());
  if ( staTrack->p() > 1.0 ) {
    mom = dirVector; 
  }

  //Get Mu state on inner muon surface
  TrajectoryStateOnSurface muTSOS = tsTransform.innerStateOnSurface(*staTrack,*theService->trackingGeometry(),&*theService->magneticField());
  
  //Get Mu state on tracker bound
  StateOnTrackerBound fromInside(&*theService->propagator(stateOnTrackerOutProp));
  //muTSOS = fromInside(muFTS);


  //Get Mu global position from the MuState (muon or tracker)
  GlobalPoint tsosPos(muTSOS.globalPosition());

  //Get the error of the global position
  GlobalError tsosErr(muTSOS.cartesianError().position());

  //Add the error to the state
  GlobalPoint tsosPosErr(tsosPos.x() + sqrt(tsosErr.cxx()),
			 tsosPos.y() + sqrt(tsosErr.cyy()),
			 tsosPos.z() + sqrt(tsosErr.czz()) );

  //Get dEta and dPhi
  float eta1 = tsosPos.eta();
  float eta2 = tsosPosErr.eta();
  float deta(fabs(eta1- eta2));
  float dphi(fabs(Geom::Phi<float>(tsosPos.phi())-Geom::Phi<float>(tsosPosErr.phi())));

  LogDebug("Analyzer") << "Region2 deta, dphi: " << deta << " " << dphi;  

  //Get Mu FTS vertex position  
  GlobalPoint vertexPos = (muFTS.position());
  GlobalError vertexErr = (muFTS.cartesianError().position());
  
  double minPt    = max(1.5,mom.perp()*0.6);
  double deltaZ   = min(15.9,3*sqrt(vertexErr.czz()));
  
  double deltaEta = 0.05;//0.05
  double deltaPhi = 0.07;//0.07

  /*    
  if ( deta > 0.05 ) { // 0.06
    deltaEta += deta/2;
  }
  if ( dphi > 0.07 ) {
    deltaPhi += 0.15;
    if ( fabs(eta1) < 1.0 && mom.perp() < 6. ) deltaPhi = dphi;
  }
  if ( fabs(eta1) < 1.25 && fabs(eta1) > 0.8 ) deltaEta = max(0.07,deltaEta);
  if ( fabs(eta1) < 1.3  && fabs(eta1) > 1.0 ) deltaPhi = max(0.3,deltaPhi);
  
  deltaEta = 1 * max(double(1. * deta),deltaEta);
  deltaPhi = 1 * max(double(1. * dphi),deltaPhi);
*/  
  
  deltaEta = deta;
  deltaPhi = dphi;

  RectangularEtaPhiTrackingRegion rectRegion(dirVector, vertexPos,
                                             minPt, 0.2,
                                             deltaZ, deltaEta, deltaPhi);
  h_region2->Fill(deltaEta,deltaPhi);
  h_region2a->Fill(deta,dphi);
  return rectRegion;

}

//
// define a region of interest within the tracker
//
RectangularEtaPhiTrackingRegion TrackingRegionAnalyzer::defineRegionOfInterest3(const reco::TrackRef& staTrack) const {

  //region 3 is FTS momentum and TSOS momentum dEta dPhi

  TrajectoryStateTransform tsTransform;
  FreeTrajectoryState muFTS = tsTransform.initialFreeState(*staTrack,&*theService->magneticField());
  
  //Get Track direction at vertex
  GlobalVector dirVector(muFTS.momentum());

  //Get track momentum
  const math::XYZVector& mo = staTrack->innerMomentum();
  GlobalVector mom(mo.x(),mo.y(),mo.z());
  if ( staTrack->p() > 1.0 ) {
    mom = dirVector; 
  }

  //-- Choices
  //-- 1) error from freeInitialState ***
  //-- 2) error from tracker bound
  //-- 3) error from inner muon surface

  //Get Mu state on inner muon surface
  //TrajectoryStateOnSurface muTSOS = tsTransform.innerStateOnSurface(*staTrack,*theService->trackingGeometry(),&*theService->magneticField());
  
  //Get Mu state on tracker bound
  //StateOnTrackerBound fromInside(&*theService->propagator(stateOnTrackerOutProp));
  //muTSOS = fromInside(muFTS);

  GlobalError  dirErr(muFTS.cartesianError().matrix().sub(4,6));
  GlobalVector dirVecErr(dirVector.x() + sqrt(dirErr.cxx()),
			 dirVector.y() + sqrt(dirErr.cyy()),
			 dirVector.z() + sqrt(dirErr.czz()));
  
  //Get dEta and dPhi
  float eta1 = dirVector.eta();
  float eta2 = dirVecErr.eta();
  float deta(fabs(eta1- eta2));
  float dphi(fabs(Geom::Phi<float>(dirVector.phi())-Geom::Phi<float>(dirVecErr.phi())));

  //LogDebug("Analyzer") << "&&&&&&& Sta" << trk->eta() << " " << trk->phi() << " " << trk->pt();
  LogDebug("Analyzer") << "Region3 deta, dphi: " << deta << " " << dphi;  
  
  GlobalPoint vertexPos = (muFTS.position());
  GlobalError vertexErr = (muFTS.cartesianError().position());
  
  double minPt    = max(1.5,mom.perp()*0.6);
  double deltaZ   = min(15.9,3*sqrt(vertexErr.czz()));
  
  double deltaEta = 0.1;//0.05
  double deltaPhi = 0.1;//0.07

   
  if ( deta > 0.05 ) {
    deltaEta += deta/2;
  }
  if ( dphi > 0.07 ) {
    deltaPhi += 0.15;
  }

  deltaPhi = min(double(0.2), deltaPhi);
  if(mom.perp() < 25.) deltaPhi = max(double(dphi),0.3);
  if(mom.perp() < 10.) deltaPhi = max(deltaPhi,0.8);
 
  deltaEta = min(double(0.2), deltaEta);
  if( mom.perp() < 6.0 ) deltaEta = 0.5;
  if( fabs(eta1) > 2.25 ) deltaEta = 0.6;
  if( fabs(eta1) > 3.0 ) deltaEta = 1.0;
  //if( fabs(eta1) > 2. && mom.perp() < 10. ) deltaEta = 1.;
  //if ( fabs(eta1) < 1.25 && fabs(eta1) > 0.8 ) deltaEta= max(0.07,deltaEta);
  if ( fabs(eta1) < 1.3  && fabs(eta1) > 1.0 ) deltaPhi = max(0.3,deltaPhi);

  //deltaEta = min(double(0.2),deltaEta);

  deltaEta = min(double(1.), 1.25 * deltaEta);
  deltaPhi = 1.2 * deltaPhi;
    
  LogDebug("Analyzer")<< "Region3 dEta, dPhi: " << deltaEta << " " << deltaPhi;
  
  //-- Choices
  //-- 4) position on inner muon surface ***
  //-- 5) position on tracker bound
  
  //Get innerMu position
  const math::XYZPoint& po = staTrack->innerPosition();
  GlobalPoint pos(po.x(),po.y(),po.z());
  
  //pos = muTSOS.globalPosition();
  
  float eta3 = pos.eta();
  float deta2(fabs(eta1- eta3));
  float dphi2(fabs(Geom::Phi<float>(dirVector.phi())-Geom::Phi<float>(pos.phi())));
  
  LogDebug("Analyzer") << "Region3 deta2, dphi2: " << deta2 << " " << dphi2;  
  
  double deltaEta2 = 0.05;
  double deltaPhi2 = 0.07;
  
  if ( deta2 > 0.05 ) { // 0.06
    deltaEta2 += deta2 / 2;
  }
  if ( dphi2 > 0.07 ) {
    deltaPhi2 += 0.15;
    if ( fabs(eta3) < 1.0 && mom.perp() < 6. ) deltaPhi2 = dphi2;
  }
  if ( fabs(eta1) < 1.25 && fabs(eta1) > 0.8 ) deltaEta2=max(0.07,deltaEta2);
  if ( fabs(eta1) < 1.3  && fabs(eta1) > 1.0 ) deltaPhi2=max(0.3,deltaPhi2);
  
  deltaEta2 = 1 * max(double(2.5 * deta2),deltaEta2);
  deltaPhi2 = 1 * max(double(3.5 * dphi2),deltaPhi2);
  
  LogDebug("Analyzer")<< "Region3 dEta2, dPhi2: " << deltaEta2 << " " << deltaPhi2;
  
  deltaEta = min(deltaEta,deltaEta2);
  deltaPhi = min(deltaPhi,deltaPhi2);
  
  RectangularEtaPhiTrackingRegion rectRegion(dirVector, vertexPos,
                                             minPt, 0.2,
                                             deltaZ, deltaEta, deltaPhi);
  h_region3->Fill(deltaEta,deltaPhi);
  h_region3a->Fill(deta,dphi);
  h_region3b->Fill(deta2,dphi2);
  return rectRegion;
  
}

//
// select tracks within the region of interest
//
reco::TrackCollection
TrackingRegionAnalyzer::chooseRegionalTrackerTracks(const RectangularEtaPhiTrackingRegion& regionOfInterest, 
                                                         const reco::TrackCollection& tkTs) const {
  
  //Get region's etaRange and phiMargin
  Range etaRange = regionOfInterest.etaRange();
  Margin phiMargin = regionOfInterest.phiMargin();

  reco::TrackCollection result;

  TrackCollection::const_iterator is;
  for ( is = tkTs.begin(); is != tkTs.end(); ++is ) {
    //check if each trackCand is in region of interest
    bool inEtaRange = etaRange.inside(is->eta());
    bool inPhiRange = (fabs(Geom::Phi<float>(is->phi()) - Geom::Phi<float>(regionOfInterest.direction().phi())) < phiMargin.right() ) ? true : false ;

    //for each trackCand in region, add trajectory and add to result
    if( inEtaRange && inPhiRange ) result.push_back(*is);
    
  }
  
  return result; 
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackingRegionAnalyzer);
