// -*- C++ -*-
//
// Package:    ChargeFilter
// Class:      ChargeFilter
// 
/**\class ChargeFilter ChargeFilter.cc UserCode/ChargeFilter/src/ChargeFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam Everett
//         Created:  Mon Oct 27 20:59:07 CDT 2008
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

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "PhysicsTools/Utilities/interface/ChargeSelector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

using namespace edm;
using namespace reco;

//
// class declaration
//

class ChargeFilter : public edm::EDFilter {
   public:
      explicit ChargeFilter(const edm::ParameterSet&);
      ~ChargeFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  int charge_;
  edm::InputTag simLabel_;
  edm::InputTag trkLabel_;
  edm::InputTag trkTPAssocLabel_;
  ChargeSelector *chargeSelector_;
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
ChargeFilter::ChargeFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  trkTPAssocLabel_ = iConfig.getParameter<InputTag>("trkMuAssocLabel");
  charge_ = iConfig.getParameter<int>("charge");
  //chargeSelector_ = new ChargeSelector(charge);
  simLabel_ = iConfig.getParameter<InputTag>("simLabel");
  trkLabel_ = iConfig.getParameter<InputTag>("trkLabel");
}


ChargeFilter::~ChargeFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ChargeFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
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


   SimToRecoCollection simToTrkColl;
   RecoToSimCollection trkToSimColl;

   Handle<SimToRecoCollection> simToTrkMuHandle;
   bool simToRecAvail = iEvent.getByLabel(trkTPAssocLabel_, simToTrkMuHandle);
   simToTrkColl = *(simToTrkMuHandle.product());
   
   Handle<RecoToSimCollection> trkMuToSimHandle;
   bool recToSimAvail = iEvent.getByLabel(trkTPAssocLabel_, trkMuToSimHandle);
   trkToSimColl = *(trkMuToSimHandle.product());

   edm::Handle<View<Track> >  trackCollection;
   iEvent.getByLabel(trkLabel_, trackCollection);

   bool goodTk = false;
   bool goodTP = false;
   bool goodMatch = false;

   ChargeSelector chargeSelector_(charge_);
   std::cout << "size " << trackCollection->size() << std::endl;
   for(View<Track>::size_type i=0; i<trackCollection->size(); ++i){
     RefToBase<Track> track(trackCollection, i);
     goodMatch = false;
     goodTk = chargeSelector_(*track);
     //     std::cout << "track i " << i << " charge " << track->charge() << " goodTk " << goodTk << std::endl;
     std::vector<std::pair<TrackingParticleRef, double> > tp;
     if(trkToSimColl.find(track) != trkToSimColl.end()){
       tp = trkToSimColl[track];
       goodTP = chargeSelector_(*tp.begin()->first);
       //       std::cout << "track i " << i << " charge " << tp.begin()->first->charge() << " goodTP " << goodTP << std::endl;
       goodMatch = goodTk & goodTP;
       //       std::cout << "goodMatch " << goodMatch;
     }
   }

   return goodTk;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ChargeFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ChargeFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ChargeFilter);
