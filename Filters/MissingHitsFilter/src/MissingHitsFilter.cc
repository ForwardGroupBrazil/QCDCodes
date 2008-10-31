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
// $Id: MissingHitsFilter.cc,v 1.1 2008/10/20 20:41:09 aeverett Exp $
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

class MissingHitsFilter : public edm::EDFilter {
   public:
      explicit MissingHitsFilter(const edm::ParameterSet&);
      ~MissingHitsFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::InputTag muonLabel_;
  int hitCut_;
  double fractionCut_;
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
  fractionCut_ = iConfig.getParameter<double>("hitFraction");
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

   /*
   for (reco::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++muon) {
     if(muon->isGlobalMuon()) {
       if(fabs(muon->eta()) > 1.2) returnVal = true;
     }
   }

   return returnVal;
   */

   
   for (reco::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++muon) {
     if(muon->isGlobalMuon()) {
       //returnVal = false;
       int hitTk = muon->track().get()->recHitsSize();
       int hitGlbTk =   muon->combinedMuon().get()->hitPattern().numberOfValidTrackerHits();
       int hitSta =  muon->standAloneMuon().get()->recHitsSize();
       int hitGlbSta = muon->combinedMuon().get()->hitPattern().numberOfValidMuonHits();
       int hitGlb =  muon->combinedMuon().get()->hitPattern().numberOfValidHits(
);

       int missingSta = hitSta-hitGlbSta;
       double fractionSta = (hitSta > 0) ? (hitSta-hitGlbSta)/hitSta : 0;
       int missingTk = hitTk-hitGlbTk;
       double fractionTk = (hitTk > 0) ? (hitTk-hitGlbTk)/hitTk : 0;
       
       if(hitCut_ > 0) {
	 if(missingSta >= hitCut_ || missingTk >= hitCut_) returnVal = true;
       }

       //if(fabs(muon->combinedMuon().get()->outerZ()) > 1100) returnVal = true;

       //else if( fractionSta >= fractionCut_ || fractionTk >= fractionCut_) returnVal = true;

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

//define this as a plug-in
DEFINE_FWK_MODULE(MissingHitsFilter);
