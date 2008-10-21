// -*- C++ -*-
//
// Package:    KaonDecay
// Class:      KaonDecay
// 
/**\class KaonDecay KaonDecay.cc Analyzer/KaonDecay/src/KaonDecay.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Chang Liu
//         Created:  Thu Sep 18 13:03:34 EDT 2008
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"

#include "TrackingTools/GeomPropagators/interface/TrackerBounds.h"
#include "TrackingTools/GeomPropagators/interface/MuonBounds.h"


//
// class declaration
//

class KaonDecay : public edm::EDFilter {
   public:
      explicit KaonDecay(const edm::ParameterSet&);
      ~KaonDecay();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  bool inTrackerBound_;
  bool inMuonBound_;
  bool inCal_;
  bool simMuon_;
};

//
// constants, enums and typedefs
//
using namespace std;
using namespace edm;

//
// static data member definitions
//

//
// constructors and destructor
//
KaonDecay::KaonDecay(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  inTrackerBound_ = iConfig.getUntrackedParameter<bool>("inTrackerBound",false);
  inMuonBound_ = iConfig.getUntrackedParameter<bool>("inMuonBound",false);
  inCal_ = iConfig.getUntrackedParameter<bool>("inCal",false);
  simMuon_ =  iConfig.getUntrackedParameter<bool>("simMuon",false);
}


KaonDecay::~KaonDecay()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
KaonDecay::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

     Handle<edm::SimTrackContainer> simTracks;
     iEvent.getByLabel("g4SimHits",simTracks);

     Handle<edm::SimVertexContainer> simVertexs;
     iEvent.getByLabel("g4SimHits",simVertexs);

     int partype = 0;

     double decayR = 0;
     double decayZ = 0;

     bool inTracker = false;
     bool inMuon = false;

     bool simMuon = false;

     for (SimTrackContainer::const_iterator p = simTracks->begin(); p != simTracks->end(); ++p){

        if (abs((*p).type()) == 13  ) {
	   simMuon = true;

           float px = p->momentum().x();
           float py = p->momentum().y();
           float pz = p->momentum().z();
           GlobalVector mom(px,py,pz);

           LogDebug("KaonDecay") << "G4 Muon: mom "<<mom<<endl;
           LogDebug("KaonDecay") << "G4 Muon: vertIndex "<<p->vertIndex()<<endl;
//           LogDebug("KaonDecay") << "G4 Muon: genpartIndex "<<p->genpartIndex()<<endl;
           if ( p->vertIndex() >= 0 && p->vertIndex() < simVertexs->size() ) {
               decayR = (*simVertexs)[p->vertIndex()].position().Rho();
               decayZ = (*simVertexs)[p->vertIndex()].position().z();
	       GlobalPoint pos((*simVertexs)[p->vertIndex()].position().x(),(*simVertexs)[p->vertIndex()].position().y(),(*simVertexs)[p->vertIndex()].position().z());
	       if(TrackerBounds::isInside(pos)) inTracker = true;
	       if(MuonBounds::isInside(pos)) inMuon = true;
               LogDebug("KaonDecay") << "G4 Vertex: pos "<<(*simVertexs)[p->vertIndex()].position() << " in tracker: " << inTracker;
               LogDebug("KaonDecay") << "G4 Vertex: parentIndex "<<(*simVertexs)[p->vertIndex()].parentIndex()<<endl;
               if ( (*simVertexs)[p->vertIndex()].parentIndex() >0 && (*simVertexs)[p->vertIndex()].parentIndex() < simTracks->size() ) {
                   partype = (*simTracks)[(*simVertexs)[p->vertIndex()].parentIndex()].type();
                   LogDebug("KaonDecay") << "G4 Vertex: parent "<<partype <<" mom "<<  (*simTracks)[(*simVertexs)[p->vertIndex()].parentIndex()].momentum()<<endl;
               }

           }
        }
     }

     if (simMuon_ && simMuon) return true;
     
     if ( partype >100 ) { //the decay-on-flight particle id, 111 or 211 or ... 
       // one can also require the decay vertex position in certain region 
       //just use decayR and decayZ
       
       LogDebug("KaonDecay") << "foundone "<<endl;
       if(inTrackerBound_ && inTracker) return true;
       else if(inMuonBound_ && inMuon) return true;
       else if(inCal_ && inMuon && !inTracker) return true;
       else return false;
       //return true;
     } else {
       return false;
     } 
}

// ------------ method called once each job just before starting event loop  ------------
void 
KaonDecay::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
KaonDecay::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(KaonDecay);
