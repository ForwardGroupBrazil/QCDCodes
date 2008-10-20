// -*- C++ -*-
//
// Package:    MissingTracker
// Class:      MissingTracker
// 
/**\class MissingTracker MissingTracker.cc UserCode/MissingTracker/src/MissingTracker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam EVERETT
//         Created:  Fri May 30 14:00:49 CEST 2008
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <TH1F.h>
#include <TH2F.h>


//
// class decleration
//

class MissingTracker : public edm::EDAnalyzer {
public:
  //  typedef TrackingRecHit::ConstRecHitContainer ConstRecHitContainer;
  
public:
  explicit MissingTracker(const edm::ParameterSet&);
  ~MissingTracker();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  //  void printHits(const TrackingRecHit::TrackingRecHitCollection& hits) const;
  // ----------member data ---------------------------
  edm::InputTag muonLabel_;
  edm::InputTag glbTrackLabel_;

  double  min, max;
  int nint;
  double  minHit, maxHit;
  int nintHit;


  DQMStore* dbe_;
  std::string dirName_;
  std::string out;


  MonitorElement * h_hitsTk, * h_hitsSta, * h_hitsGlb;
  MonitorElement * h_hitsGlbTk, * h_hitsGlbSta;
  MonitorElement * h_lostHitsTk, * h_lostHitsSta;
  MonitorElement * h_hitsTk_eta, * h_hitsSta_eta, * h_hitsGlb_eta;
  MonitorElement * h_hitsGlbTk_eta, * h_hitsGlbSta_eta;
  MonitorElement * h_lostHitsTk_eta, * h_lostHitsSta_eta;

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
MissingTracker::MissingTracker(const edm::ParameterSet& iConfig) :
   min(iConfig.getParameter<double>("min")),
   max(iConfig.getParameter<double>("max")),
   nint(iConfig.getParameter<int>("nint")),
   minHit(iConfig.getParameter<double>("minHit")),
   maxHit(iConfig.getParameter<double>("maxHit")),
   nintHit(iConfig.getParameter<int>("nintHit"))
{
  //now do what ever initialization is needed
  muonLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("muLabel");
  glbTrackLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("glbLabel");

  dbe_ = edm::Service<DQMStore>().operator->();
  out = iConfig.getParameter<std::string>("out");
  dirName_ = iConfig.getParameter<std::string>("dirName");
  
}


MissingTracker::~MissingTracker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
MissingTracker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

   LogTrace("MissingTracker")<<"Event number: " << iEvent.id().event();

   //LogDebug("MissingTracker");  
   Handle<reco::MuonCollection> muons;
   iEvent.getByLabel(muonLabel_,muons);
   //   const reco::TrackCollection muons = *(muonHandle.product());

   reco::TrackRef ttrk;
   reco::TrackRef strk;
   reco::TrackRef gtrk;
   
   int hitTk = 0;
   int hitGlbTk = 0;
   int hitSta = 0;
   int hitGlbSta = 0;
   int hitGlb = 0;

   for (reco::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++muon) {
     if(muon->isGlobalMuon()) {
       ttrk = muon->track();
       strk = muon->standAloneMuon();
       gtrk = muon->combinedMuon();
       //
       
       if( ttrk.isAvailable() ) { 	
	 hitTk = ttrk.get()->hitPattern().numberOfValidTrackerHits();
       } else  edm::LogWarning("MissingTracker")<<"No ttrk";     
       
       if( strk.isAvailable() ) { 
	 hitSta = strk.get()->hitPattern().numberOfValidMuonHits();
       } else  edm::LogWarning("MissingTracker")<<"No strk";     
       
       if( gtrk.isAvailable() ) { 
	 hitGlbTk = gtrk.get()->hitPattern().numberOfValidTrackerHits();
	 hitGlbSta = gtrk.get()->hitPattern().numberOfValidMuonHits();
	 hitGlb = gtrk.get()->hitPattern().numberOfValidHits();
       } else  edm::LogWarning("MissingTracker")<<"No gtrk";     

       ////-----
       
       h_hitsTk->Fill(hitTk);
       h_hitsGlbTk->Fill(hitGlbTk);
       h_hitsSta->Fill(hitSta);
       h_hitsGlbSta->Fill(hitGlbSta);
       h_hitsGlb->Fill(hitGlb);
       
       h_lostHitsTk->Fill(hitTk-hitGlbTk);
       h_lostHitsSta->Fill(hitSta-hitGlbSta);
       
       h_hitsTk_eta->Fill(muon->eta(),hitTk);
       h_hitsGlbTk_eta->Fill(muon->eta(),hitGlbTk);
       h_hitsSta_eta->Fill(muon->eta(),hitSta);
       h_hitsGlbSta_eta->Fill(muon->eta(),hitGlbSta);
       h_hitsGlb_eta->Fill(muon->eta(),hitGlb);
       
       h_lostHitsTk_eta->Fill(muon->eta(),hitTk-hitGlbTk);
       h_lostHitsSta_eta->Fill(muon->eta(),hitSta-hitGlbSta);
       
       
       LogTrace("MissingTracker")<<"Slava    Tk  " << hitTk;
       LogTrace("MissingTracker")<<"Slava GlbTk  " << hitGlbTk;
       LogTrace("MissingTracker")<<"Slava deltaTk  " << hitTk - hitGlbTk;
       LogTrace("MissingTracker")<<"Slava    Mu  " << hitSta;
       LogTrace("MissingTracker")<<"Slava GlbMu  " << hitGlbSta;
       LogTrace("MissingTracker")<<"Slava deltaMu  " << hitSta - hitGlbSta;
     }

     //}
     // ADAM
     
     
     if(gtrk.isAvailable()) LogTrace("MissingTracker") << "Number of Glb Hits: " << gtrk->recHitsSize() << std::endl;
     //     ConstRecHitContainer muonRecHits;
     for (trackingRecHit_iterator hit = muon->combinedMuon()->recHitsBegin();  hit != muon->combinedMuon()->recHitsEnd();  ++hit) {
       //       muonRecHits.push_back(*hit);
       //       const GlobalPoint& pos = (*hit)->globalPosition();
       int idNum =  (*hit)->geographicalId();
       if((*hit)->isValid()) LogTrace("MissingTracker")
			       //	 << "r = " << sqrt(pos.x() * pos.x() + pos.y() * pos.y())
			       //	 << "  z = " << pos.z()
       	 << "  dimension = " << (*hit)->dimension()
	 << "  " << idNum
	 << "  " << (*hit)->geographicalId().det()
	 << "  " << (*hit)->geographicalId().subdetId();
     }    
     
     if(strk.isAvailable()) LogTrace("MissingTracker") << "Number of Sta Hits: " << strk->recHitsSize() << std::endl;
     for (trackingRecHit_iterator hit = muon->standAloneMuon()->recHitsBegin();  hit != muon->standAloneMuon()->recHitsEnd();  ++hit) {
       //       muonRecHits.push_back(*hit);
       //       const GlobalPoint& pos = (*hit)->globalPosition();
       int idNum =  (*hit)->geographicalId();
       if((*hit)->isValid()) LogTrace("MissingTracker")
	 //	 << "r = " << sqrt(pos.x() * pos.x() + pos.y() * pos.y())
	 //	 << "  z = " << pos.z()
       	 << "  dimension = " << (*hit)->dimension()
	 << "  " << idNum
	 << "  " << (*hit)->geographicalId().det()
	 << "  " << (*hit)->geographicalId().subdetId();
     }
     //     printHits(muonRecHits);

     if(ttrk.isAvailable()) LogTrace("MissingTracker") << "Number of Tk Hits: " << ttrk->recHitsSize() << std::endl;
     for (trackingRecHit_iterator hit = muon->track()->recHitsBegin();  hit != muon->track()->recHitsEnd();  ++hit) {
       //       const GlobalPoint& pos = (*hit)->globalPosition();
       int idNum =  (*hit)->geographicalId();
       if((*hit)->isValid()) LogTrace("MissingTracker") 
	 //	 << "r = " << sqrt(pos.x() * pos.x() + pos.y() * pos.y())
	 //	 << "  z = " << pos.z()
       	 << "  dimension = " << (*hit)->dimension()
       	 << "  " << idNum
       	 << "  " << (*hit)->geographicalId().det()
       	 << "  " << (*hit)->geographicalId().subdetId();
       else LogTrace("MissingTracker") << "invalid hit";
     }
     //LogDebug("MissingTracker");  
     // end ADAM
   }

   ///-----
   /*

   //LogDebug("MissingTracker");  
   //
   // Jim's code
   //
   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(glbTrackLabel_, tracks);
   //LogDebug("MissingTracker");  
   for (reco::TrackCollection::const_iterator track = tracks->begin();  track != tracks->end();  ++track) {   //LogDebug("MissingTracker");  
     //LogTrace("MissingTracker") << "track eta = " << track->eta();
     bool trackerhits = false;
     for (trackingRecHit_iterator hit = track->recHitsBegin();  hit != track->recHitsEnd();  ++hit) {
       if ((*hit)->geographicalId().det() == DetId::Tracker) {
	 trackerhits = true;
       }
       
       //       std::cout << "    hit id = " << (*hit)->geographicalId() << " pos = " << (*hit)->localPosition() << std::endl;
     }
     //LogDebug("MissingTracker");  
     if (!trackerhits) {
       edm::LogWarning("MissingTracker") << "This one!!!" << std::endl;
          //LogDebug("MissingTracker");  
       for (trackingRecHit_iterator hit = track->recHitsBegin();  hit != track->recHitsEnd();  ++hit) {   //LogDebug("MissingTracker");  
	 int idNum = (*hit)->geographicalId();
         LogTrace("MissingTracker") << "    hit id = " << idNum ;
       }
     }
        //LogDebug("MissingTracker");  
   }
*/
   //
   // end Jim's code
   //
   //LogDebug("MissingTracker");  
}


// ------------ method called once each job just before starting event loop  ------------
void 
MissingTracker::beginJob(const edm::EventSetup&)
{

  dbe_->cd();
  //InputTag algo = theTracksLabel;
  std::string dirName=dirName_;
  //if (algo.process()!="")
  //  dirName+=algo.process()+"_";
  //if(algo.label()!="")
  //  dirName+=algo.label()+"_";
  //if(algo.instance()!="")
  //  dirName+=algo.instance()+"";      
  //if (dirName.find("Tracks")<dirName.length()){
  //  dirName.replace(dirName.find("Tracks"),6,"");
  //}
  std::replace(dirName.begin(), dirName.end(), ':', '_');
  dbe_->setCurrentFolder(dirName.c_str());
  

  h_hitsTk= dbe_->book1D("hitsTk", "number of hits per track", nintHit,minHit,maxHit ) ;
  h_hitsSta= dbe_->book1D("hitsSta", "number of hits per track", nintHit,minHit,maxHit ) ;
  h_hitsGlb= dbe_->book1D("hitsGlb", "number of hits per track", nintHit,minHit,maxHit ) ;
  h_hitsGlbTk= dbe_->book1D("hitsGlbTk", "number of hits per track", nintHit,minHit,maxHit ) ;
  h_hitsGlbSta= dbe_->book1D("hitsGlbSta", "number of hits per track", nintHit,minHit,maxHit ) ;
  
  h_lostHitsTk= dbe_->book1D("lostHitsTk", "number of hits per track", nintHit,minHit,maxHit ) ;
  h_lostHitsSta= dbe_->book1D("lostHitsSta", "number of hits per track", nintHit,minHit,maxHit ) ;
  //h_lostHitsGlb= dbe_->book1D("lostHitsGlb", "number of hits per track", nintHit,minHit,maxHit ) ;
  
  h_hitsTk_eta= dbe_->book2D("hitsTk_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  h_hitsSta_eta= dbe_->book2D("hitsSta_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  h_hitsGlbTk_eta= dbe_->book2D("hitsGlbTk_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  h_hitsGlbSta_eta= dbe_->book2D("hitsGlbSta_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  h_hitsGlb_eta= dbe_->book2D("hitsGlb_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  h_lostHitsTk_eta= dbe_->book2D("lostHitsTk_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  h_lostHitsSta_eta= dbe_->book2D("lostHitsSta_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  //h_lostHitsGlb_eta= dbe_->book2D("lostHitsGlb_eta", "number of hits per track",nint,min,max, nintHit,minHit,maxHit ) ;
  
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MissingTracker::endJob() {
  if ( out.size() != 0 && dbe_ ) dbe_->save(out);
  
}
/*
  void MissingTracker::printHits(const TrackingRecHit::TrackingRecHitCollection& hits) const {
  
  LogTrace(theCategory) << "Used RecHits: " << hits.size();
  for (ConstRecHitContainer::const_iterator ir = hits.begin(); ir != hits.end(); ir++ ) {
  if ( !(*ir)->isValid() ) {
      LogTrace(theCategory) << "invalid RecHit";
      continue; 
    }
    
    const GlobalPoint& pos = (*ir)->globalPosition();
    
    LogTrace(theCategory) 
      << "r = " << sqrt(pos.x() * pos.x() + pos.y() * pos.y())
      << "  z = " << pos.z()
      << "  dimension = " << (*ir)->dimension()
      << "  " << (*ir)->det()->geographicalId().det()
      << "  " << (*ir)->det()->subDetector();
  }

}
*/
//define this as a plug-in
DEFINE_FWK_MODULE(MissingTracker);
