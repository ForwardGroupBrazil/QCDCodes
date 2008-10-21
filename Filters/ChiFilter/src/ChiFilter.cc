// -*- C++ -*-
//
// Package:    ChiFilter
// Class:      ChiFilter
// 
/**\class ChiFilter ChiFilter.cc UserCode/ChiFilter/src/ChiFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam A Everett
//         Created:  Wed Sep 24 13:00:00 EDT 2008
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
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//
// class declaration
//

class ChiFilter : public edm::EDFilter {
   public:
      explicit ChiFilter(const edm::ParameterSet&);
      ~ChiFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::InputTag muonLabel_;
  int chi2Cut_;
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
ChiFilter::ChiFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  muonLabel_ = iConfig.getParameter<edm::InputTag>("muLabel");
  chi2Cut_ =  iConfig.getParameter<int>("chi2Cut");
}


ChiFilter::~ChiFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ChiFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
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

   Handle<reco::MuonCollection> muons;
   iEvent.getByLabel(muonLabel_,muons);
   
   bool returnVal = false;

   for (reco::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++muon) {
     if(muon->isGlobalMuon()) {
       double chiVal = muon->combinedMuon().get()->normalizedChi2();
       if (chi2Cut_ > 0) { 
	 if(chiVal < fabs(chi2Cut_) ) returnVal = true;
       } else {
	 if(chiVal > fabs(chi2Cut_) ) returnVal = true;
       }
     }
   }

   return returnVal;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ChiFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ChiFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ChiFilter);
