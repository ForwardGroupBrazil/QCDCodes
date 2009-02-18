// -*- C++ -*-
//
// Package:    MissingHitsFilter
// Class:      MissingHitsFilter
// 
/**\class MissingHitsFilter MissingHitsFilter.cc UserCode/MissingHitsFilter/src/MissingHitsFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam A Everett
//         Created:  Wed Sep 10 15:08:42 EDT 2008
// $Id: MissingHitsFilter.cc,v 1.3 2008/11/07 16:03:34 aeverett Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

//
// class declaration
//

class MissingHitsFilter : public edm::EDFilter {
   public:
      explicit MissingHitsFilter(const edm::ParameterSet&);
      ~MissingHitsFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual int countMuonHits(const reco::Track& track) const;
      virtual int countTrackerHits(const reco::Track& track) const;
      // ----------member data ---------------------------
  edm::InputTag muonLabel_;
  int hitCut_;
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
MissingHitsFilter::MissingHitsFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  muonLabel_ = iConfig.getParameter<edm::InputTag>("muLabel");
  hitCut_ =  iConfig.getParameter<int>("hitCut");
}


MissingHitsFilter::~MissingHitsFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MissingHitsFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::MuonCollection> muons;
   iEvent.getByLabel(muonLabel_,muons);
   
   bool returnVal = false;

   for (reco::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++muon) {
     int hitTk =0;
     int hitSta =0;
     int hitGlbTk = 0;
     int hitGlbMu = 0;
     int hitGlb = 0;
     if ( muon->isTrackerMuon() ) {
       hitTk = countTrackerHits(*muon->track());
     }
     if ( muon->isStandAloneMuon() ) {
       hitSta =  countMuonHits(*muon->standAloneMuon());
     }
     if ( muon->isGlobalMuon() ) {
       hitGlbTk = countTrackerHits(*muon->combinedMuon());       
       hitGlbMu = countMuonHits(*muon->combinedMuon());
       hitGlb =  muon->combinedMuon().get()->recHitsSize();
       
       int missingSta = hitSta-hitGlbMu;
       int missingTk = hitTk-hitGlbTk;
       
       if(hitCut_ > 0) {
	 if( (missingSta >= hitCut_) || (missingTk >= hitCut_) ) {
	   returnVal = true;
	 }
       } else if (hitCut_ < 0) { 
	 if( (missingSta <= hitCut_) || (missingTk <= hitCut_) ) {
	   returnVal = true; 
	 }
       }
       
     }
   }

   return returnVal;
   
}

// ------------ method called once each job just before starting event loop  ------------
void 
MissingHitsFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MissingHitsFilter::endJob() {
}

int
MissingHitsFilter::countMuonHits(const reco::Track& track) const {
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
MissingHitsFilter::countTrackerHits(const reco::Track& track) const {
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

//define this as a plug-in
DEFINE_FWK_MODULE(MissingHitsFilter);
