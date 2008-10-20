// -*- C++ -*-
//
// Package:    MuonTrackSelector
// Class:      MuonTrackSelector
// 
/**\class MuonTrackSelector MuonTrackSelector.cc Validation/MuonTrackSelector/src/MuonTrackSelector.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Arun Luthra
//         Created:  Mon Jul 14 01:13:27 CEST 2008
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "RecoMuon/MuonIdentification/interface/IdGlobalFunctions.h"

#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"


class MuonTrackSelector : public edm::EDProducer {
public:
  explicit MuonTrackSelector(const edm::ParameterSet&);
  ~MuonTrackSelector();
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::InputTag muonsTag;
  edm::InputTag selectionTag;
  const edm::ParameterSet parset_;      
};

MuonTrackSelector::MuonTrackSelector(const edm::ParameterSet& parset) :
  muonsTag(    parset.getParameter< edm::InputTag >("muonsTag")),
  selectionTag(parset.getParameter< edm::InputTag >("selectionTag")),
  parset_(parset)
{

  
  //produces<reco::TrackCollection>();//.setBranchAlias("isGoodInnerTracks");
  //  produces<reco::TrackCollection>("outerTracks").setBranchAlias("isGoodOuterTracks");
  //  produces<reco::TrackCollection>("globalTracks").setBranchAlias("isGoodGlobalTracks");

  //aaa
  std::string alias( parset.getParameter<std::string>( "@module_label" ) );  
  produces<reco::TrackCollection>("TrackerOnly").setBranchAlias( alias + "TrackerOnlyTracks" );
  produces<reco::TrackExtraCollection>("TrackerOnly").setBranchAlias( alias + "TrackerOnlyExtras" );
  produces<TrackingRecHitCollection>("TrackerOnly").setBranchAlias( alias + "TrackerOnlyHits" );
 
}


MuonTrackSelector::~MuonTrackSelector()
{
}

void
MuonTrackSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace reco;

  edm::Handle<reco::MuonCollection> muonCollectionH;
  
  edm::LogVerbatim("MuonTrackSelector") << std::endl << "muonsTag = " << muonsTag.label() << std::endl;
  edm::LogVerbatim("MuonTrackSelector") << "selectionTag = " << selectionTag.label() << std::endl;
  
  iEvent.getByLabel(muonsTag,muonCollectionH);
  if(muonCollectionH->size()==0) {
    edm::LogVerbatim("MuonTrackSelector") << "\nNo muons in this event\n";
  }
  
  //std::auto_ptr<reco::TrackCollection> goodTracks(new reco::TrackCollection);
  
  //adam
  
  std::auto_ptr<reco::TrackCollection> selTracks_(new reco::TrackCollection);
  std::auto_ptr<reco::TrackExtraCollection> selTracksExtras_( new reco::TrackExtraCollection );
  std::auto_ptr<TrackingRecHitCollection> selTracksHits_( new TrackingRecHitCollection );
  
  TrackingRecHitRefProd rHits = iEvent.getRefBeforePut<TrackingRecHitCollection>("TrackerOnly");
  reco::TrackExtraRefProd rTrackExtras = iEvent.getRefBeforePut<TrackExtraCollection>("TrackerOnly");
  reco::TrackRefProd rTracks = iEvent.getRefBeforePut<TrackCollection>("TrackerOnly");    
  
  //end adam
  
  size_t hidx = 0, idx = 0;  
  for(reco::MuonCollection::const_iterator muon = muonCollectionH->begin(); muon != muonCollectionH->end(); ++muon) {
    bool isGoodSelection = false;
    
    if (selectionTag.label() == "AllTrackerMuons") {
      isGoodSelection = muon->isGood(reco::Muon::AllTrackerMuons);
      if(isGoodSelection) {
	//goodTracks->push_back( *muon->track() );
	
	// adam
	TrackRef trkRef = muon->track();
	if(trkRef.isNonnull()){
	  
	  selTracks_->push_back(Track( *trkRef) );
	  
	  Track & trk= selTracks_->back();
	  
	  selTracksExtras_->push_back( TrackExtra( trk.outerPosition(), trk.outerMomentum(), trk.outerOk(),
						   trk.innerPosition(), trk.innerMomentum(), trk.innerOk(),
						   trk.outerStateCovariance(), trk.outerDetId(),
						   trk.innerStateCovariance(), trk.innerDetId(),
						   trk.seedDirection() ) );
	  
	  TrackExtra & tx = selTracksExtras_->back();
	  
	  for( trackingRecHit_iterator hit = trk.recHitsBegin(); hit != trk.recHitsEnd(); ++ hit ) {
	    selTracksHits_->push_back( (*hit)->clone() );
	    tx.add( TrackingRecHitRef( rHits, hidx ++ ) );
	  }
	  
	  trk.setExtra( TrackExtraRef( rTrackExtras, idx ++ ) );
	  
	}
	//end adam
      }
    }
  }
  //iEvent.put(goodTracks);

  //adam
     iEvent.put( selTracks_ , "TrackerOnly");
     iEvent.put( selTracksExtras_ , "TrackerOnly");
     iEvent.put( selTracksHits_ ,"TrackerOnly");
     
}

void 
MuonTrackSelector::beginJob(const edm::EventSetup&)
{
}

void 
MuonTrackSelector::endJob() {
}

DEFINE_FWK_MODULE(MuonTrackSelector);
