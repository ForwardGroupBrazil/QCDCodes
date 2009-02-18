// -*- C++ -*-
//
// Package:    MissingTracker
// Class:      MissingTracker
// 
/**\class MissingTracker MissingTracker.cc UserCode/MissingTracker/src/MissingTracker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam EVERETT
//         Created:  Fri May 30 14:00:49 CEST 2008
// $Id: MissingTracker.cc,v 1.1 2008/10/20 20:44:13 aeverett Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"



#include <TH1F.h>
#include <TH2F.h>


//
// class decleration
//

class MissingTracker : public edm::EDAnalyzer {
public:
  //  typedef TrackingRecHit::ConstRecHitContainer ConstRecHitContainer;
  
public:
  explicit MissingTracker(const edm::ParameterSet&);
  ~MissingTracker();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void setEvent(const edm::Event& event);
  virtual int countTrackerHits(const reco::Track& track) const;
  virtual int countMuonHits(const reco::Track& track) const;
  virtual TransientTrackingRecHit::ConstRecHitContainer getTransientRecHits(const reco::Track& track) const;
  virtual void printHits(const TransientTrackingRecHit::ConstRecHitContainer& hits) const;
  // ----------member data ---------------------------
  edm::InputTag muonLabel_;
  edm::InputTag glbTrackLabel_;

  MuonServiceProxy * theService;
  std::string theTrackerRecHitBuilderName;
  edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
  
  std::string theMuonRecHitBuilderName;
  edm::ESHandle<TransientTrackingRecHitBuilder> theMuonRecHitBuilder;

  std::string theTrackerPropagatorName;
  
  unsigned long long theCacheId_TRH;
  bool theRPCInTheFit;

  double  min, max;
  int nint;
  double  minHit, maxHit;
  int nintHit;


  DQMStore* dbe_;
  std::string dirName_;
  std::string out;


  MonitorElement * h_hitsTk, * h_hitsSta, * h_hitsGlb;
  MonitorElement * h_hitsGlbTk, * h_hitsGlbSta;
  MonitorElement * h_lostHitsTk, * h_lostHitsSta;
  MonitorElement * h_hitsTk_eta, * h_hitsSta_eta, * h_hitsGlb_eta;
  MonitorElement * h_hitsGlbTk_eta, * h_hitsGlbSta_eta;
  MonitorElement * h_lostHitsTk_eta, * h_lostHitsSta_eta;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MissingTracker::MissingTracker(const edm::ParameterSet& iConfig) :
   min(iConfig.getParameter<double>("min")),
   max(iConfig.getParameter<double>("max")),
   nint(iConfig.getParameter<int>("nint")),
   minHit(iConfig.getParameter<double>("minHit")),
   maxHit(iConfig.getParameter<double>("maxHit")),
   nintHit(iConfig.getParameter<int>("nintHit"))
{
  //now do what ever initialization is needed
  muonLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("muLabel");
  glbTrackLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("glbLabel");

  dbe_ = edm::Service<DQMStore>().operator->();
  out = iConfig.getParameter<std::string>("out");
  dirName_ = iConfig.getParameter<std::string>("dirName");
  
  // the service parameters
  edm::ParameterSet serviceParameters 
    = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);

  theRPCInTheFit = iConfig.getParameter<bool>("RefitRPCHits");

  theTrackerPropagatorName = iConfig.getParameter<std::string>("TrackerPropagator");

  theTrackerRecHitBuilderName = iConfig.getParameter<std::string>("TrackerRecHitBuilder");
  theMuonRecHitBuilderName = iConfig.getParameter<std::string>("MuonRecHitBuilder");  

  theCacheId_TRH = 0;

}


MissingTracker::~MissingTracker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  if ( theService ) delete theService;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
MissingTracker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   setEvent(iEvent);

   LogTrace("MissingTracker")<<"Event number: " << iEvent.id().event();

   //LogDebug("MissingTracker");  
   Handle<reco::MuonCollection> muons;
   iEvent.getByLabel(muonLabel_,muons);
   //   const reco::TrackCollection muons = *(muonHandle.product());


   
   for (reco::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++muon) {
     reco::TrackRef ttrk;
     reco::TrackRef strk;
     reco::TrackRef gtrk;
     
     if(muon->isGlobalMuon()) {
       int hitTk = 0;
       int hitGlbTk = 0;
       int hitSta = 0;
       int hitGlbSta = 0;
       int hitGlb = 0;
       
       ttrk = muon->track();
       strk = muon->standAloneMuon();
       gtrk = muon->combinedMuon();
       
       if( ttrk.isAvailable() ) { 	
	 hitTk = countTrackerHits(*ttrk);
       } else  edm::LogWarning("MissingTracker")<<"No ttrk";     
       
       if( strk.isAvailable() ) { 
	 hitSta = countMuonHits(*strk);
       } else  edm::LogWarning("MissingTracker")<<"No strk";     
       
       if( gtrk.isAvailable() ) { 
	 hitGlbTk = countTrackerHits(*gtrk);
	 hitGlbSta = countMuonHits(*gtrk);
	 hitGlb = gtrk->recHitsSize();
	 //gtrk.get()->hitPattern().print();
       } else  edm::LogWarning("MissingTracker")<<"No gtrk";     

       h_hitsTk->Fill(hitTk);
       h_hitsGlbTk->Fill(hitGlbTk);
       h_hitsSta->Fill(hitSta);
       h_hitsGlbSta->Fill(hitGlbSta);
       h_hitsGlb->Fill(hitGlb);
       
       h_lostHitsTk->Fill(hitTk-hitGlbTk);
       h_lostHitsSta->Fill(hitSta-hitGlbSta);
       
       h_hitsTk_eta->Fill(muon->eta(),hitTk);
       h_hitsGlbTk_eta->Fill(muon->eta(),hitGlbTk);
       h_hitsSta_eta->Fill(muon->eta(),hitSta);
       h_hitsGlbSta_eta->Fill(muon->eta(),hitGlbSta);
       h_hitsGlb_eta->Fill(muon->eta(),hitGlb);
       
       h_lostHitsTk_eta->Fill(muon->eta(),hitTk-hitGlbTk);
       h_lostHitsSta_eta->Fill(muon->eta(),hitSta-hitGlbSta);
              
       LogTrace("MissingTracker")<<"Slava      Tk  " << hitTk;
       LogTrace("MissingTracker")<<"Slava   GlbTk  " << hitGlbTk;
       LogTrace("MissingTracker")<<"Slava deltaTk  " << hitTk - hitGlbTk;
       LogTrace("MissingTracker")<<"Slava      Mu  " << hitSta;
       LogTrace("MissingTracker")<<"Slava   GlbMu  " << hitGlbSta;
       LogTrace("MissingTracker")<<"Slava deltaMu  " << hitSta - hitGlbSta;
       LogTrace("MissingTracker")<<"Slava     Glb  " << hitGlb;
 
       if(ttrk.isAvailable()) {
	 LogTrace("MissingTracker") << "Number of Tk Hits: " << ttrk->recHitsSize() << std::endl;
	 printHits(getTransientRecHits(*ttrk));
       }
       if(strk.isAvailable()) {
	 LogTrace("MissingTracker") << "Number of Sta Hits: " << strk->recHitsSize() << std::endl;
	 printHits(getTransientRecHits(*strk));     
       }
       if(gtrk.isAvailable()) {
	 LogTrace("MissingTracker") << "Number of Glb Hits: " << gtrk->recHitsSize() << std::endl;
	 printHits(getTransientRecHits(*gtrk));     
       }
     } //end isGlobalMuon


   }// end loop over muon collection

}


// ------------ method called once each job just before starting event loop  ------------
void 
MissingTracker::beginJob(const edm::EventSetup& eventSetup)
{
  if ( theService ) theService->update(eventSetup);
  dbe_->cd();
  //InputTag algo = theTracksLabel;
  std::string dirName=dirName_;
  //if (algo.process()!="")
  //  dirName+=algo.process()+"_";
  //if(algo.label()!="")
  //  dirName+=algo.label()+"_";
  //if(algo.instance()!="")
  //  dirName+=algo.instance()+"";      
  //if (dirName.find("Tracks")<dirName.length()){
  //  dirName.replace(dirName.find("Tracks"),6,"");
  //}
  std::replace(dirName.begin(), dirName.end(), ':', '_');
  dbe_->setCurrentFolder(dirName.c_str());
  

  h_hitsTk= dbe_->book1D("hitsTk", "number of hits per track", nintHit,minHit,maxHit ) ;
  h_hitsSta= dbe_->book1D("hitsSta", "number of hits per track", nintHit,minHit,maxHit ) ;
  h_hitsGlb= dbe_->book1D("hitsGlb", "number of hits per track", nintHit,minHit,maxHit ) ;
  h_hitsGlbTk= dbe_->book1D("hitsGlbTk", "number of hits per track", nintHit,minHit,maxHit ) ;
  h_hitsGlbSta= dbe_->book1D("hitsGlbSta", "number of hits per track", nintHit,minHit,maxHit ) ;
  
  h_lostHitsTk= dbe_->book1D("lostHitsTk", "number of hits per track", nintHit,minHit,maxHit ) ;
  h_lostHitsSta= dbe_->book1D("lostHitsSta", "number of hits per track", nintHit,minHit,maxHit ) ;
  //h_lostHitsGlb= dbe_->book1D("lostHitsGlb", "number of hits per track", nintHit,minHit,maxHit ) ;
  
  h_hitsTk_eta= dbe_->book2D("hitsTk_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  h_hitsSta_eta= dbe_->book2D("hitsSta_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  h_hitsGlbTk_eta= dbe_->book2D("hitsGlbTk_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  h_hitsGlbSta_eta= dbe_->book2D("hitsGlbSta_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  h_hitsGlb_eta= dbe_->book2D("hitsGlb_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  h_lostHitsTk_eta= dbe_->book2D("lostHitsTk_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  h_lostHitsSta_eta= dbe_->book2D("lostHitsSta_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  //h_lostHitsGlb_eta= dbe_->book2D("lostHitsGlb_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MissingTracker::endJob() {
  if ( out.size() != 0 && dbe_ ) dbe_->save(out);
  
}

void MissingTracker::printHits(const TransientTrackingRecHit::ConstRecHitContainer& hits) const {
  
  LogTrace("MissingTracker") << "Used RecHits: " << hits.size();
  for (TransientTrackingRecHit::ConstRecHitContainer::const_iterator ir = hits.begin(); ir != hits.end(); ir++ ) {
    if ( !(*ir)->isValid() ) {
      LogTrace("MissingTracker") << "invalid RecHit";
      continue; 
    }
    
    const GlobalPoint& pos = (*ir)->globalPosition();
    
    LogTrace("MissingTracker") 
      << "r = " << sqrt(pos.x() * pos.x() + pos.y() * pos.y())
      << "  z = " << pos.z()
      << "  R = " << sqrt(pos.x() * pos.x() + pos.y() * pos.y() + pos.z() * pos.z())
      << " type = " << (*ir)->type()
      << "  dimension = " << (*ir)->dimension()
      << "  " << (*ir)->det()->geographicalId().det()
      << "  " << (*ir)->det()->subDetector();
  } 
} 


int
MissingTracker::countMuonHits(const reco::Track& track) const {
  TransientTrackingRecHit::ConstRecHitContainer result;
  
  int count = 0;

  for (trackingRecHit_iterator hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
    if((*hit)->isValid()) {
      DetId recoid = (*hit)->geographicalId();
      if ( recoid.det() == DetId::Muon ) count++;
    }
  }
  return count;
}

int
MissingTracker::countTrackerHits(const reco::Track& track) const {
  TransientTrackingRecHit::ConstRecHitContainer result;
  
  int count = 0;

  for (trackingRecHit_iterator hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
    if((*hit)->isValid()) {
      DetId recoid = (*hit)->geographicalId();
      if ( recoid.det() == DetId::Tracker ) count++;
    }
  }
  return count;
}

TransientTrackingRecHit::ConstRecHitContainer
MissingTracker::getTransientRecHits(const reco::Track& track) const {

  TransientTrackingRecHit::ConstRecHitContainer result;
  
  TrajectoryStateTransform tsTrans;

  TrajectoryStateOnSurface currTsos = tsTrans.innerStateOnSurface(track, *theService->trackingGeometry(), &*theService->magneticField());

  for (trackingRecHit_iterator hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
    if((*hit)->isValid()) {
      DetId recoid = (*hit)->geographicalId();
      if ( recoid.det() == DetId::Tracker ) {
	TransientTrackingRecHit::RecHitPointer ttrhit = theTrackerRecHitBuilder->build(&**hit);
	TrajectoryStateOnSurface predTsos =  theService->propagator(theTrackerPropagatorName)->propagate(currTsos, theService->trackingGeometry()->idToDet(recoid)->surface());
	//LogTrace(theCategory)<<"predtsos "<<predTsos.isValid();
	if ( predTsos.isValid() ) currTsos = predTsos;
	TransientTrackingRecHit::RecHitPointer preciseHit = ttrhit->clone(predTsos);
	result.push_back(preciseHit);
      } else if ( recoid.det() == DetId::Muon ) {
	if ( (*hit)->geographicalId().subdetId() == 3 && !theRPCInTheFit) {
	  //LogDebug(theCategory) << "RPC Rec Hit discarded"; 
	  continue;
	}
	result.push_back(theMuonRecHitBuilder->build(&**hit));
      }
    }
  }
  
  return result;

  }

void MissingTracker::setEvent(const edm::Event& event) {

  unsigned long long newCacheId_TRH = theService->eventSetup().get<TransientRecHitRecord>().cacheIdentifier();
  if ( newCacheId_TRH != theCacheId_TRH ) {
    theCacheId_TRH = newCacheId_TRH;
    theService->eventSetup().get<TransientRecHitRecord>().get(theTrackerRecHitBuilderName,theTrackerRecHitBuilder);
    theService->eventSetup().get<TransientRecHitRecord>().get(theMuonRecHitBuilderName,theMuonRecHitBuilder);

  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(MissingTracker);
