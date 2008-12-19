/**  \class GlobalTruncator
 * 
 *   TeV muon reconstructor:
 *
 *
 *   $Date: 2008/12/18 18:04:16 $
 *   $Revision: 1.2 $
 *
 *   \author  Piotr Traczyk (SINS Warsaw)
 *   \author  Adam Everett (Purdue University)
 */

// Framework
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "UserCode/GlobalTruncator/src/GlobalTruncator.h"

// TrackFinder and specific GLB Trajectory Builder
#include "RecoMuon/GlobalTrackFinder/interface/GlobalMuonTrajectoryBuilder.h"
#include "RecoMuon/TrackingTools/interface/MuonTrackFinder.h"
#include "RecoMuon/TrackingTools/interface/MuonTrackLoader.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/GlobalTrackingTools/interface/GlobalMuonRefitter.h"

// Input and output collection
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"

using namespace edm;
using namespace std;

//
// constructor with config
//
GlobalTruncator::GlobalTruncator(const ParameterSet& parameterSet) {

  LogDebug("Muon|RecoMuon|GlobalTruncator") << "constructor called" << endl;

  // GLB Muon Collection Label
  theGLBCollectionLabel = parameterSet.getParameter<InputTag>("MuonCollectionLabel");

  // service parameters
  ParameterSet serviceParameters = parameterSet.getParameter<ParameterSet>("ServiceParameters");

  // the services
  theService = new MuonServiceProxy(serviceParameters);
  
  // TrackRefitter parameters
  ParameterSet refitterParameters = parameterSet.getParameter<ParameterSet>("RefitterParameters");
  theRefitter = new GlobalTruncRefitter(refitterParameters, theService);

  // TrackLoader parameters
  ParameterSet trackLoaderParameters = parameterSet.getParameter<ParameterSet>("TrackLoaderParameters");
  theTrackLoader = new MuonTrackLoader(trackLoaderParameters,theService);

  theRefits = parameterSet.getParameter< std::vector<std::string> >("Refits");
  theRefitIndex = parameterSet.getParameter< std::vector<int> >("RefitIndex");
  theRefitSubIndex = parameterSet.getParameter< std::vector<int> >("RefitSubIndex");

  for(unsigned int ww=0;ww<theRefits.size();ww++){
    LogDebug("Muon|RecoMuon|GlobalTruncator") << "Refit " << theRefits[ww];
    produces<reco::TrackCollection>(theRefits[ww]);
    produces<TrackingRecHitCollection>(theRefits[ww]);
    produces<reco::TrackExtraCollection>(theRefits[ww]);
    produces<vector<Trajectory> >(theRefits[ww]) ;
    produces<TrajTrackAssociationCollection>(theRefits[ww]);
    produces<reco::TrackToTrackMap>(theRefits[ww]);
  }
}


//
// destructor
//
GlobalTruncator::~GlobalTruncator() {

  LogTrace("Muon|RecoMuon|GlobalTruncator") << "destructor called" << endl;
  if (theService) delete theService;
  if (theRefitter) delete theRefitter;

}


//
// reconstruct muons
//
void GlobalTruncator::produce(Event& event, const EventSetup& eventSetup) {

  const string metname = "Muon|RecoMuon|GlobalTruncator";  
  LogTrace(metname)<< endl << endl;
  LogTrace(metname)<< "TeV Muon Reconstruction started" << endl;  

  // Update the services
  theService->update(eventSetup);

  theRefitter->setEvent(event);

  theRefitter->setServices(theService->eventSetup());

  // Take the GLB muon container(s)
  Handle<reco::TrackCollection> glbMuons;
  event.getByLabel(theGLBCollectionLabel,glbMuons);

  Handle<vector<Trajectory> > glbMuonsTraj;

  LogTrace(metname)<< "Taking " << glbMuons->size() << " Global Muons "<<theGLBCollectionLabel<<endl;

  vector<MuonTrajectoryBuilder::TrackCand> glbTrackCands;

  event.getByLabel(theGLBCollectionLabel.label(), glbMuonsTraj);
    
  const reco::TrackCollection *glbTracks = glbMuons.product();
  
  for(unsigned int ww=0;ww<theRefits.size();ww++) {
    LogDebug(metname)<<"TeVRefit for Refit: " <<theRefitIndex[ww];
    std::vector<std::pair<Trajectory*,reco::TrackRef> > miniMap;
    vector<Trajectory*> trajectories;
    reco::TrackRef::key_type trackIndex = 0;
    for (reco::TrackCollection::const_iterator track = glbTracks->begin(); track!=glbTracks->end(); track++ , ++trackIndex) {
      reco::TrackRef glbRef(glbMuons,trackIndex);
      
      vector<Trajectory> refitted=theRefitter->refit(*track,theRefitIndex[ww],theRefitSubIndex[ww]);

      if (refitted.size()>0) {
        Trajectory *refit = new Trajectory(refitted.front());
	LogDebug(metname)<<"TeVTrackLoader for Refit: " <<theRefits[ww];
	trajectories.push_back(refit);
	std::pair<Trajectory*,reco::TrackRef> thisPair(refit,glbRef);
	miniMap.push_back(thisPair);
      }
    }
    theTrackLoader->loadTracks(trajectories,event,miniMap,theRefits[ww]);
  }
    
  LogTrace(metname) << "Done." << endl;    

}