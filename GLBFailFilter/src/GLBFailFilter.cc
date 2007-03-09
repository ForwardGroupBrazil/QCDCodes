// -*- C++ -*-
//
// Package:    GLBFailFilter
// Class:      GLBFailFilter
// 
/**\class GLBFailFilter GLBFailFilter.cc UserCode/GLBFailFilter/src/GLBFailFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam Everett
//         Created:  Fri Mar  9 21:26:27 CET 2007
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

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"



//
// class declaration
//

class GLBFailFilter : public edm::EDFilter {
   public:
      explicit GLBFailFilter(const edm::ParameterSet&);
      ~GLBFailFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::InputTag TKtrackTags_; 
  edm::InputTag STAtrackTags_; 
  edm::InputTag MuonTags_; 
  edm::InputTag SIMtrackTags_; 

  double PtCut_;

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
GLBFailFilter::GLBFailFilter(const edm::ParameterSet& iConfig)
  :
  TKtrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("TKtracks")),
  STAtrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("STAtracks")),
  MuonTags_(iConfig.getUntrackedParameter<edm::InputTag>("Muons")),
  SIMtrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("SIMtracks")),
  PtCut_(iConfig.getParameter<double>("PtCut"))
{
  //now do what ever initialization is needed
  
}


GLBFailFilter::~GLBFailFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
GLBFailFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;

  bool result = true;
  
  //using namespace edm;
  using reco::TrackCollection;
  using reco::MuonCollection;
  
  Handle<edm::SimTrackContainer> SIMTrackCollection;
  iEvent.getByLabel( SIMtrackTags_,SIMTrackCollection);
  const SimTrackContainer simTC = *(SIMTrackCollection.product());
  
  Handle<reco::TrackCollection> TKTrackCollection;  
  iEvent.getByLabel( TKtrackTags_, TKTrackCollection);
  const reco::TrackCollection tkTC = *(TKTrackCollection.product());

  Handle<reco::TrackCollection> STATrackCollection;    
  iEvent.getByLabel( STAtrackTags_, STATrackCollection);
  const reco::TrackCollection staTC = *(STATrackCollection.product());
  
  Handle<reco::MuonCollection> MuCollection;  
  iEvent.getByLabel( MuonTags_, MuCollection);
  const reco::MuonCollection muonC = *(MuCollection.product());
  
  int nGoodSTA = 0;
  TrackCollection::const_iterator staTrack1;
  for(staTrack1 = staTC.begin(); staTrack1 != staTC.end(); ++staTrack1){
    if ( ! ((*staTrack1).pt() < PtCut_ 
	    || (*staTrack1).innerMomentum().Rho() < PtCut_ 
	    || (*staTrack1).innerMomentum().R() < 2.5 ) ) 
      nGoodSTA++;
  }
  
  bool goodSTA = (nGoodSTA > 0) ? true : false ;
  
  if(muonC.size() == 0 && goodSTA && tkTC.size() > 0) result = false;



   return result;
}

// ------------ method called once each job just before starting event loop  ------------
void 
GLBFailFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GLBFailFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(GLBFailFilter);
