// -*- C++ -*-
//
// Package:    PtFilter
// Class:      PtFilter
// 
/**\class PtFilter PtFilter.cc UserCode/PtFilter/src/PtFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam A Everett
//         Created:  Wed Sep 10 15:08:42 EDT 2008
// $Id: PtFilter.cc,v 1.3 2008/11/07 16:03:34 aeverett Exp $
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

//
// class declaration
//

class PtFilter : public edm::EDFilter {
   public:
      explicit PtFilter(const edm::ParameterSet&);
      ~PtFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::InputTag muonLabel_;
  double maxPtCut_;
  double minPtCut_;
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
PtFilter::PtFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  muonLabel_ = iConfig.getParameter<edm::InputTag>("muLabel");
  maxPtCut_ =  iConfig.getParameter<double>("MaxPtCut");
  minPtCut_ =  iConfig.getParameter<double>("MinPtCut");

}


PtFilter::~PtFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
PtFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

   using namespace edm;

   Handle<reco::MuonCollection> muons;
   iEvent.getByLabel(muonLabel_,muons);
   
   bool returnVal = false;
   
   for (reco::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++muon) {
     if( muon->isGlobalMuon() ) {
       if(muon->pt() > minPtCut_ && muon->pt() < maxPtCut_) returnVal = true;
     }
   }
   
   return returnVal;
      
}

// ------------ method called once each job just before starting event loop  ------------
void 
PtFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PtFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PtFilter);
