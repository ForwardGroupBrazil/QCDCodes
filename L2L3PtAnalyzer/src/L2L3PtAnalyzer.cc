// -*- C++ -*-
//
// Package:    L2L3PtAnalyzer
// Class:      L2L3PtAnalyzer
// 
/**\class L2L3PtAnalyzer L2L3PtAnalyzer.cc UserCode/L2L3PtAnalyzer/src/L2L3PtAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam A Everett
//         Created:  Tue Jan 20 14:12:47 EST 2009
// $Id: L2L3PtAnalyzer.cc,v 1.1.2.2 2009/01/20 21:14:18 aeverett Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"

//
// class decleration
//

class L2L3PtAnalyzer : public edm::EDAnalyzer {
   public:
      explicit L2L3PtAnalyzer(const edm::ParameterSet&);
      ~L2L3PtAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  reco::TrackCollection getSeededTkCollection(const reco::TrackRef&, const edm::View<reco::Track>& ) ;
      // ----------member data ---------------------------
  edm::InputTag L2Label_;
  edm::InputTag L3TkLabel_;
  edm::InputTag TPLabel_;
  edm::InputTag trkMuAssocLabel_;

  TrackAssociatorBase* trkMuAssociator_;

};

//
// constants, enums and typedefs
//
using namespace std;
using namespace edm;
using namespace reco;

//
// static data member definitions
//

//
// constructors and destructor
//
L2L3PtAnalyzer::L2L3PtAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  L2Label_ = iConfig.getParameter<InputTag>("L2Label");
  L3TkLabel_ = iConfig.getParameter<InputTag>("L3TkLabel");
  TPLabel_ = iConfig.getParameter<InputTag>("TPLabel");
  trkMuAssocLabel_ = iConfig.getParameter<InputTag>("trkMuAssocLabel");

}


L2L3PtAnalyzer::~L2L3PtAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
L2L3PtAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle<TrackingParticleCollection>  TPCollectionH ;
   iEvent.getByLabel(TPLabel_,TPCollectionH);
   const TrackingParticleCollection tPC = *(TPCollectionH.product());
   
   edm::Handle<TrackCollection>  L2Collection;
   iEvent.getByLabel(L2Label_, L2Collection);

   LogDebug("L2L3PtAnalyzer")<<"L2Collection size "<< L2Collection->size();

   edm::Handle<View<Track> >  L3TkCollection;
   iEvent.getByLabel(L3TkLabel_, L3TkCollection);
   View<Track> L3TkColl = *(L3TkCollection.product());

   LogDebug("L2L3PtAnalyzer")<<"L3TkCollection size "<< L3TkCollection->size();

   reco::RecoToSimCollection recSimColl;
   reco::SimToRecoCollection simRecColl;
   
   LogTrace("L2L3PtAnalyzer") << "Calling associateRecoToSim method" << "\n";
   recSimColl=trkMuAssociator_->associateRecoToSim(L3TkCollection,
						   TPCollectionH,
						   &iEvent);
   
   LogTrace("L2L3PtAnalyzer") << "Calling associateSimToReco method" << "\n";
   simRecColl=trkMuAssociator_->associateSimToReco(L3TkCollection,
						   TPCollectionH, 
						   &iEvent);
   
   for(TrackCollection::size_type i=0; i<L2Collection->size(); ++i){
     reco::TrackRef staTrack(L2Collection,i);     
     TrackCollection seededTkCollection = getSeededTkCollection(staTrack,L3TkColl);     
     LogDebug("L2L3PtAnalyzer")<<"SeededTkCollection size "<< seededTkCollection.size();
     
   }
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
L2L3PtAnalyzer::beginJob(const edm::EventSetup& eventSetup)
{
  ESHandle<TrackAssociatorBase> trkMuAssocHandle;
  eventSetup.get<TrackAssociatorRecord>().get(trkMuAssocLabel_.label(), trkMuAssocHandle);
  trkMuAssociator_ = const_cast<TrackAssociatorBase*>(trkMuAssocHandle.product());
  
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L2L3PtAnalyzer::endJob() {
}

TrackCollection
L2L3PtAnalyzer::getSeededTkCollection(const reco::TrackRef& staTrack, const View<Track>& L3TkCollection ) {
  TrackCollection tkTrackCands;

  for(View<Track>::const_iterator iTk=L3TkCollection.begin(); iTk != L3TkCollection.end() ; ++iTk) {

    //reco::TrackRef tkTrack(L3TkCollection,j);
    edm::Ref<L3MuonTrajectorySeedCollection> l3seedRef = iTk->seedRef().castTo<edm::Ref<L3MuonTrajectorySeedCollection> >() ;
    reco::TrackRef staTrack_2 = l3seedRef->l2Track();
    if(staTrack_2 == staTrack ) { 
      tkTrackCands.push_back(*iTk);
    }
  }
  return tkTrackCands;
}

//define this as a plug-in
DEFINE_FWK_MODULE(L2L3PtAnalyzer);
