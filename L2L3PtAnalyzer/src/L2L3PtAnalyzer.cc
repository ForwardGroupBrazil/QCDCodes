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
// $Id: L2L3PtAnalyzer.cc,v 1.1.2.4 2009/01/20 22:25:28 aeverett Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h" 

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

#include "PhysicsTools/Utilities/interface/PtComparator.h"
#include <TH1.h>
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
  std::vector<edm::RefToBase<reco::Track> > getSeededTkCollection(const reco::TrackRef&, const edm::Handle<edm::View<reco::Track> >& ) ;

  struct MyPtComparator {
    bool operator()(const edm::RefToBase<reco::Track>& a,
		    const edm::RefToBase<reco::Track>& b) const {
      return a->pt() > b->pt(); 
    }
  };

      // ----------member data ---------------------------
  edm::InputTag L2Label_;
  edm::InputTag L3TkLabel_;
  edm::InputTag TPLabel_;
  edm::InputTag trkMuAssocLabel_;

  TrackAssociatorBase* trkMuAssociator_;

  TH1 * h_pdg_1;
  TH1 * h_pdg_2;
  TH1 * h_pdg_3;
  TH1 * h_pdg_4;

  TH1 * h_Pt_sta;

  TH1 * h_Pt_1;
  TH1 * h_Pt_2;
  TH1 * h_Pt_3;
  TH1 * h_Pt_4;

  TH1 * h_deltaPt_1;
  TH1 * h_deltaPt_2;
  TH1 * h_deltaPt_3;
  TH1 * h_deltaPt_4;

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
     vector<RefToBase<Track> > seededTkCollection = getSeededTkCollection(staTrack,L3TkCollection);     
     LogDebug("L2L3PtAnalyzer")<<"SeededTkCollection size "<< seededTkCollection.size();
     
     h_Pt_sta->Fill(fabs(staTrack->pt()));

     //sort by Pt
     //GreaterByPt<RefToBase<Track> > compTracks;
     stable_sort(seededTkCollection.begin(),seededTkCollection.end(),MyPtComparator());
     
     //for each seededTk, get the associated TP
     int iParticle = 0; 
     for(vector<RefToBase<Track> >::const_iterator track = seededTkCollection.begin(); track != seededTkCollection.end(); ++track){
       iParticle++;

       if(iParticle==1) h_Pt_1->Fill(fabs((*track)->pt()));
       if(iParticle==2) h_Pt_2->Fill(fabs((*track)->pt()));
       if(iParticle==3) h_Pt_3->Fill(fabs((*track)->pt()));
       if(iParticle==4) h_Pt_4->Fill(fabs((*track)->pt()));

       if(iParticle==1) h_deltaPt_1->Fill(fabs((*track)->pt()- staTrack->pt()));
       if(iParticle==2) h_deltaPt_2->Fill(fabs((*track)->pt()- staTrack->pt()));
       if(iParticle==3) h_deltaPt_3->Fill(fabs((*track)->pt()- staTrack->pt()));
       if(iParticle==4) h_deltaPt_4->Fill(fabs((*track)->pt()- staTrack->pt()));

       std::vector<std::pair<TrackingParticleRef, double> > tp;
       if(recSimColl.find(*track) != recSimColl.end()){
	 tp = recSimColl[*track];
	 if (tp.size()!=0) {
	   edm::LogVerbatim("L2L3PtAnalyzer") << "reco::Track #" 
					      << iParticle << " with pt=" 
					      << (*track)->pt() 
					      << " associated with quality:" 
					      << tp.begin()->second 
					      <<" to tpgId " 
					      << tp.begin()->first->pdgId() <<"\n";
	   if(iParticle==1) h_pdg_1->Fill(fabs(tp.begin()->first->pdgId()));
	   if(iParticle==2) h_pdg_2->Fill(fabs(tp.begin()->first->pdgId()));
	   if(iParticle==3) h_pdg_3->Fill(fabs(tp.begin()->first->pdgId()));
	   if(iParticle==4) h_pdg_4->Fill(fabs(tp.begin()->first->pdgId()));
	 }
       } else {
	 edm::LogVerbatim("L2L3PtAnalyzer") << "reco::Track #" 
					    << iParticle << " with pt=" 
					    << (*track)->pt()
					    << " NOT associated to any TrackingParticle" << "\n";		  
	   if(iParticle==1) h_pdg_1->Fill(333);
	   if(iParticle==2) h_pdg_2->Fill(333);
	   if(iParticle==3) h_pdg_3->Fill(333);
	   if(iParticle==4) h_pdg_4->Fill(333);
       }       
     }     
   }
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
L2L3PtAnalyzer::beginJob(const edm::EventSetup& eventSetup)
{
  ESHandle<TrackAssociatorBase> trkMuAssocHandle;
  eventSetup.get<TrackAssociatorRecord>().get(trkMuAssocLabel_.label(), trkMuAssocHandle);
  trkMuAssociator_ = const_cast<TrackAssociatorBase*>(trkMuAssocHandle.product());
  
  edm::Service<TFileService> fs;
  h_pdg_1 = fs->make<TH1F>("h_pdg_1","PDGId of Track 1",501,-0.5,500.5);
  h_pdg_2 = fs->make<TH1F>("h_pdg_2","PDGId of Track 1",501,-0.5,500.5);
  h_pdg_3 = fs->make<TH1F>("h_pdg_3","PDGId of Track 1",501,-0.5,500.5);
  h_pdg_4 = fs->make<TH1F>("h_pdg_4","PDGId of Track 1",501,-0.5,500.5);

  h_Pt_sta=fs->make<TH1F>("h_Pt_sta","p_{T}^{#mu}",100,0.,1000.);

  h_Pt_1=fs->make<TH1F>("h_Pt_1","p_{T}^{1}",100,0.,1000.);
  h_Pt_2=fs->make<TH1F>("h_Pt_2","p_{T}^{2}",100,0.,1000.);
  h_Pt_3=fs->make<TH1F>("h_Pt_3","p_{T}^{3}",100,0.,1000.);
  h_Pt_4=fs->make<TH1F>("h_Pt_4","p_{T}^{4}",100,0.,1000.);

  h_deltaPt_1=fs->make<TH1F>("h_deltaPt_1","#Delta(p_{T}^{#mu,1})",100,0.,1000.);
  h_deltaPt_2=fs->make<TH1F>("h_deltaPt_2","#Delta(p_{T}^{#mu,2})",100,0.,1000.);
  h_deltaPt_3=fs->make<TH1F>("h_deltaPt_3","#Delta(p_{T}^{#mu,3})",100,0.,1000.);
  h_deltaPt_4=fs->make<TH1F>("h_deltaPt_4","#Delta(p_{T}^{#mu,4})",100,0.,1000.);
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L2L3PtAnalyzer::endJob() {
}

vector<RefToBase<Track> >
L2L3PtAnalyzer::getSeededTkCollection(const reco::TrackRef& staTrack, const Handle<View<Track> >& L3TkCollection ) {
  TrackCollection tkTrackCands;
  vector<RefToBase<Track> > tkTrackRefs;
  
  for(View<Track>::size_type i=0; i<L3TkCollection->size(); ++i){
    RefToBase<Track> track(L3TkCollection, i);
    edm::Ref<L3MuonTrajectorySeedCollection> l3seedRef = track->seedRef().castTo<edm::Ref<L3MuonTrajectorySeedCollection> >() ;
    reco::TrackRef staTrack_2 = l3seedRef->l2Track();
    if(staTrack_2 == staTrack ) { 
      tkTrackRefs.push_back(track);
    }
  }
  return tkTrackRefs;
}

//define this as a plug-in
DEFINE_FWK_MODULE(L2L3PtAnalyzer);
