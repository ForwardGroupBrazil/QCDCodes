// -*- C++ -*-
//
// Package:    GlobalMatchingAnalyser
// Class:      GlobalMatchingAnalyser
// 
/**\class GlobalMatchingAnalyser GlobalMatchingAnalyser.cc FastAnalysis/GlobalMatchingAnalyser/src/GlobalMatchingAnalyser.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam Everett
//         Created:  Fri Dec 18 12:47:08 CST 2009
// $Id: GlobalMatchingAnalyser.cc,v 1.14 2010/11/03 20:59:11 aeverett Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include <DataFormats/TrackReco/interface/Track.h>
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/GlobalTrackingTools/interface/MuonTrackingRegionBuilder.h"
#include "RecoMuon/GlobalTrackingTools/interface/GlobalMuonTrackMatcher.h"
#include "RecoTracker/TkTrackingRegions/interface/RectangularEtaPhiTrackingRegion.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"


#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrackFitters/interface/RecHitLessByDet.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoMuon/GlobalTrackingTools/interface/GlobalMuonRefitter.h"
#include "RecoMuon/TrackingTools/interface/MuonCandidate.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"


#include "TFile.h"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

//
// class decleration
//

class GlobalMatchingAnalyser : public edm::EDAnalyzer {
public:
  explicit GlobalMatchingAnalyser(const edm::ParameterSet&);
  ~GlobalMatchingAnalyser();

  typedef std::pair<const Trajectory*, reco::TrackRef> TrackCand;
  typedef TransientTrackingRecHit::ConstRecHitContainer ConstRecHitContainer;

  typedef std::vector<Trajectory> TC;
  typedef TC::const_iterator TI;
  typedef MuonCandidate::CandidateContainer CandidateContainer;

protected:
  struct ComparatorInOut {
    
    bool operator()(const TransientTrackingRecHit::ConstRecHitPointer& a,
		    const TransientTrackingRecHit::ConstRecHitPointer& b) const{ 
      bool barrel_a = ( a->det()->subDetector() == GeomDetEnumerators::DT ||
			a->det()->subDetector() == GeomDetEnumerators::RPCBarrel );
      
      bool barrel_b = ( b->det()->subDetector() == GeomDetEnumerators::DT ||
			b->det()->subDetector() == GeomDetEnumerators::RPCBarrel );
      
      if ( barrel_a && barrel_b ) return  a->det()->surface().position().perp() < b->det()->surface().position().perp();
      
      else if ( !barrel_a && !barrel_b ) return  fabs(a->globalPosition().z()) < fabs(b->globalPosition().z());
      else if ( barrel_a && !barrel_b  ) return true;
      else if ( !barrel_a && barrel_b  ) return false;
      //shouldn;t really get here in any case (there's some sense to throw here )
      return false;
    }
  };
  enum RefitDirection{inToOut,outToIn,undetermined};

  /// check order of RechIts on a trajectory
  RefitDirection checkRecHitsOrdering(const ConstRecHitContainer&) const;

    

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  RectangularEtaPhiTrackingRegion defineRegionOfInterest(const reco::TrackRef&) const;
  
  std::vector<TrackCand> 
  chooseRegionalTrackerTracksFixed(const reco::TrackRef& staCand,
				   const std::vector<TrackCand>& tkTs,
				   int iSta);
  
  std::vector<TrackCand> 
  chooseRegionalTrackerTracks(const reco::TrackRef& staCand,
			      const std::vector<TrackCand>& tkTs,
			      int iSta);

  ////////

  CandidateContainer
  build(const reco::TrackRef staCand, const reco::TrackRef tkCand) const;
  
  TransientTrackingRecHit::ConstRecHitContainer
  getTransientRecHits(const reco::Track& track) const;

  double 
  trackProbability(const Trajectory& track) const;
  
  // ----------member data ---------------------------
  MuonTrackingRegionBuilder* theRegionBuilder;

  MuonServiceProxy* theService;

  GlobalMuonTrackMatcher* theTrackMatcher;

  std::map<std::string, std::map<std::string, TH1*> > testPlots;

  //TGraphErrors *tg1;
  TGraphErrors *glb_combined, *glb_inner, *glb_outer;
  TGraphErrors *sta_muon, *sta_muonsurf;
  TGraphErrors *tk_muon;

  TGraphErrors *region_dynamic, *region_fixed;

  TGraphErrors *pos_tkCand, *pos_tkCandFixed, *pos_tkAll, *pos_tkAllsurf;
  TGraphErrors *pos_selectedTkCand, *pos_selectedTkCandFixed;

  TGraphErrors *surface_error1;
  TGraphErrors *surface_error2;

  TH1F *h_distance_muHit, *h_Distance_muHit, *h_chi2_muHit, *h_loc_chi2_muHit, *h_loc_chi2_2_muHit, *h_deltaR_muHit, *h_loc_chi2_3_muHit;
  TH1F *h_distance_tkHit, *h_Distance_tkHit, *h_chi2_tkHit, *h_loc_chi2_tkHit, *h_loc_chi2_2_tkHit, *h_deltaR_tkHit, *h_loc_chi2_3_tkHit;
  TH1F *h_distance_tksurf, *h_Distance_tksurf, *h_chi2_tksurf, *h_loc_chi2_tksurf, *h_loc_chi2_2_tksurf, *h_deltaR_tksurf, *h_loc_chi2_3_tksurf;

  TH1F *h_sta_pt, *h_sta_p, *h_sta_rho, *h_sta_R, *h_sta_cut, *h_sta_etaFlip, *h_sta_etaFlip1, *h_muon_steps, *h_muon_refit, *h_muon_refit_chi2;

  TH1F *h_sta_updated, *h_muon_type;

  TH1F *h_sta_rho_outer, *h_sta_R_outer;

  TH2F *h_distance_loc_chi2_muHit;
  TH2F *h_distance_loc_chi2_tkHit;
  TH2F *h_distance_loc_chi2_tksurf;

  TH2F *h_distance_loc_chi2_2_muHit;
  TH2F *h_distance_loc_chi2_2_tkHit;
  TH2F *h_distance_loc_chi2_2_tksurf;

  TH2F *h_distance_loc_chi2_3_muHit;
  TH2F *h_distance_loc_chi2_3_tkHit;
  TH2F *h_distance_loc_chi2_3_tksurf;

  edm::InputTag theTrackLabel;
  edm::InputTag theMuonLabel;
  edm::InputTag classif_;
  edm::InputTag theSTACollectionLabel;

  int useAll;

  std::string theTrackerRecHitBuilderName;
  edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
  
  std::string theMuonRecHitBuilderName;
  edm::ESHandle<TransientTrackingRecHitBuilder> theMuonRecHitBuilder;
  
  unsigned long long theCacheId_TRH;
  std::string theTrackerPropagatorName;
  bool theRPCInTheFit;

  GlobalMuonRefitter* theGlbRefitter;
  int   theMuonHitsOption;
  TrackTransformer* theTrackTransformer;

  bool fixForFlip;

};

//
// constants, enums and typedefs
//
typedef std::pair<const Trajectory*, reco::TrackRef> TrackCand;

//
// static data member definitions
//

using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//
GlobalMatchingAnalyser::GlobalMatchingAnalyser(const edm::ParameterSet& iConfig) : theRegionBuilder(0), theService(0), theTrackMatcher(0)

{
  //using namespace edm;
  
  //now do what ever initialization is needed

  theTrackLabel = iConfig.getParameter<edm::InputTag>("trackLabel");
  theMuonLabel = iConfig.getParameter<edm::InputTag>("muonLabel");
  theSTACollectionLabel = iConfig.getParameter<edm::InputTag>("staLabel");

  useAll = iConfig.getParameter<int>("useAll");

  fixForFlip = iConfig.getUntrackedParameter<bool>("fixForFlip",false);

  classif_ = iConfig.getParameter<edm::InputTag>("classification");

  // the services
  ParameterSet serviceParameters = iConfig.getParameter<ParameterSet>("ServiceParameters");  
  theService = new MuonServiceProxy(serviceParameters);

  // the region builder
  ParameterSet regionBuilderPSet = iConfig.getParameter<ParameterSet>("MuonTrackingRegionBuilder");  
  theRegionBuilder = new MuonTrackingRegionBuilder(regionBuilderPSet,theService);

  ParameterSet trackMatcherPSet = iConfig.getParameter<ParameterSet>("GlobalMuonTrackMatcher");
  theTrackMatcher = new GlobalMuonTrackMatcher(trackMatcherPSet,theService);

  theTrackerRecHitBuilderName = iConfig.getParameter<string>("TrackerRecHitBuilder");
  theMuonRecHitBuilderName = iConfig.getParameter<string>("MuonRecHitBuilder");  
  theCacheId_TRH = 0;

  theTrackerPropagatorName = iConfig.getParameter<string>("TrackerPropagator");
  theRPCInTheFit = iConfig.getParameter<bool>("RefitRPCHits");

  // TrackRefitter parameters
  ParameterSet refitterParameters = iConfig.getParameter<ParameterSet>("GlbRefitterParameters");
  theGlbRefitter = new GlobalMuonRefitter(refitterParameters, theService);

  theMuonHitsOption = refitterParameters.getParameter<int>("MuonHitsOption");

  ParameterSet trackTransformerPSet = iConfig.getParameter<ParameterSet>("TrackTransformer");

  theTrackTransformer = new TrackTransformer(trackTransformerPSet);

}


GlobalMatchingAnalyser::~GlobalMatchingAnalyser()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  if(theRegionBuilder) delete theRegionBuilder;
  if(theService) delete theService;
  if(theTrackMatcher) delete theTrackMatcher;
  if (theTrackTransformer) delete theTrackTransformer;
  if (theGlbRefitter) delete theGlbRefitter;


}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GlobalMatchingAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //using namespace edm;
  //using namespace reco;
  
  theService->update(iSetup);
  
  theRegionBuilder->setEvent(iEvent);

  theTrackTransformer->setServices(theService->eventSetup());

  theGlbRefitter->setEvent(iEvent);
  theGlbRefitter->setServices(theService->eventSetup());

  unsigned long long newCacheId_TRH = theService->eventSetup().get<TransientRecHitRecord>().cacheIdentifier();
  if ( newCacheId_TRH != theCacheId_TRH ) {
    LogDebug("MatchAnalyzer") << "TransientRecHitRecord changed!";
    theCacheId_TRH = newCacheId_TRH;
    theService->eventSetup().get<TransientRecHitRecord>().get(theTrackerRecHitBuilderName,theTrackerRecHitBuilder);
    theService->eventSetup().get<TransientRecHitRecord>().get(theMuonRecHitBuilderName,theMuonRecHitBuilder);
  }
  
  
  LogDebug("MatchAnalyzer") << "********************" << "Run " << iEvent.id().run() << " Event " << iEvent.id().event() ;

  /////////////////////////////
  //
  // Get the collections
  //
  /////////////////////////////
  
  // Get Muons
  Handle<View<Muon> > muonHandle;
  iEvent.getByLabel(theMuonLabel, muonHandle);
  View<Muon> muonColl = *(muonHandle.product());
  
  // Get Muon Tracks
  Handle<View<Track> > trkHandle;
  iEvent.getByLabel(theTrackLabel, trkHandle);
  View<Track> trkColl = *(trkHandle.product());
  
  edm::Handle<reco::TrackCollection> allTrackerTracks;
  iEvent.getByLabel(theTrackLabel,allTrackerTracks);

  /*  
  Handle<ValueMap<int> > classif; iEvent.getByLabel(classif_, classif);

  for (size_t i = 0, n = muonColl.size(); i < n; ++i) {
    edm::RefToBase<reco::Muon> muRef = muonColl.refAt(i);
    int origin = (*classif)[muRef];
    LogDebug("MatchAnalyzer");
    if (origin < 0) {
      LogTrace("MatchAnalyzer") << "This is a ghost!" << std::endl;
    } else {
      switch (origin) {   
      case 0: LogTrace("MatchAnalyzer") << "Unmatched ()" << std::endl; break;
      case 1: LogTrace("MatchAnalyzer") << "Fake" << std::endl; break;
      case 2: LogTrace("MatchAnalyzer") << "Light flavour or decay" << std::endl; break;          
      case 3: LogTrace("MatchAnalyzer") << "Heavy flavour or tau" << std::endl; break;
      case 4: LogTrace("MatchAnalyzer") << "Primary muon" << std::endl; break;
      }
    }
  }
  */

  /////////////////////////////
  //
  // Make a collection for all tracker tracks
  //
  /////////////////////////////
  vector<TrackCand> tkTrackCands;
  
  for ( unsigned int position = 0; position != allTrackerTracks->size(); ++position ) {
    reco::TrackRef tkTrackRef(allTrackerTracks,position);
    TrackCand tkCand = TrackCand((Trajectory*)(0),tkTrackRef);
    tkTrackCands.push_back(tkCand); 
  }

  /////////////////////////////
  //
  // Set the point for all tracker tracks
  //
  /////////////////////////////
  
  int iTkAll = 0;
  for(vector<TrackCand>::const_iterator iTk = tkTrackCands.begin();
      iTk != tkTrackCands.end(); ++iTk) {
    //TrackRef tkRef(tkPreCandColl,position);
    //position++;
    pos_tkAll->SetPoint(iTkAll,iTk->second->eta(),iTk->second->phi());
    pos_tkAll->SetPointError(iTkAll,iTk->second->etaError(),iTk->second->phiError());
    iTkAll++;
    //tkCandColl.push_back(TrackCand((Trajectory*)(0),tkRef));
  }
  
  // // Test crap
  //tg1->SetPoint(0,1,1);
  //tg1->SetPoint(1,2,1);
  //tg1->SetPointError(0,1,.5);
  //tg1->SetPointError(1,0.25,1.);

  int nMu = muonColl.size();

  /////////////////////////////
  //
  // Loop over all muons
  //
  /////////////////////////////

  int iMu = 0;
  int iSta = 0;
  int iTkMu = 0;
  int iGlbMu = 0;
  int surfaceOffset = 0;
  int iTkFixed = 0;
  int iTkDynamic = 0;
  int iSelTkFixed = 0;
  int iSelTkDynamic = 0;
  bool etaFlipBool = true;
  bool propBool = true;
  for(View<Muon>::const_iterator iMuon = muonColl.begin();
      iMuon != muonColl.end(); ++iMuon, iMu++) {
    etaFlipBool = true;
    propBool = true;
    h_muon_steps->Fill(0);
    LogTrace("MatchAnalyzer") << "*****" 
			      << endl << "Muon " << iMu+1 << " of " 
			      << nMu << endl;

    LogTrace("MatchAnalyzer") << "     isXXX? " 
			      << iMuon->isGlobalMuon() << " " 
			      << iMuon->isStandAloneMuon() << " " 
			      << iMuon->isTrackerMuon() << endl;

    const reco::TrackRef glbTrack = ( iMuon->isGlobalMuon()) ? 
      iMuon->combinedMuon() : reco::TrackRef();

    const reco::TrackRef tkTrack = ( iMuon->isTrackerMuon() ) ? 
      iMuon->innerTrack() : TrackRef();
    
    const reco::TrackRef staTrack = ( iMuon->isStandAloneMuon() ) ? 
      iMuon->outerTrack() : TrackRef();

    h_muon_type->Fill(0);
    if(iMuon->isGlobalMuon()) h_muon_type->Fill(1);
    if(iMuon->isTrackerMuon()) h_muon_type->Fill(2);
    if(iMuon->isStandAloneMuon()) h_muon_type->Fill(3);
    reco::TrackRef staTrackForRegion = staTrack;
    //////////////////////////////
    //
    // Special (temporary) selection
    //
    //////////////////////////////
    //bool good = true; //(!iMuon->isGlobalMuon() && iMuon->isTrackerMuon() && iMuon->isStandAloneMuon());
    //if( !good ) continue;
    //LogTrace("MatchAnalyzer") <<"     Passes selection";
    h_muon_steps->Fill(1);
    if(glbTrack.isAvailable())  LogTrace("MatchAnalyzer") <<"                    Glb pT " << glbTrack->pt();
    if(tkTrack.isAvailable())   LogTrace("MatchAnalyzer") <<"                    Tk  pT " << tkTrack->pt();
    if(staTrack.isAvailable())  LogTrace("MatchAnalyzer") <<"                    Sta pT " << staTrack->pt();

    bool etaFliped = false;

    if(staTrack.isAvailable()) {
      h_muon_steps->Fill(2);
      h_sta_pt->Fill(staTrack->pt());
      h_sta_p->Fill(staTrack->p());
      h_sta_rho->Fill(staTrack->innerMomentum().rho());
      h_sta_R->Fill(staTrack->innerMomentum().R());
      h_sta_rho_outer->Fill(staTrack->outerMomentum().rho());
      h_sta_R_outer->Fill(staTrack->outerMomentum().R());
      h_sta_cut->Fill(0);
      int cut = (staTrack->pt() < 1. || staTrack->innerMomentum().rho() < 1. || staTrack->innerMomentum().R() < 2.5) ? -1 : 1;
      h_sta_cut->Fill(cut);
      int cut2 = (staTrack->pt() < 1. && staTrack->innerMomentum().rho() < 1. && staTrack->innerMomentum().R() < 2.5) ? -2 : 2 ;
      h_sta_cut->Fill(cut2);
      
      if(!staTrack.isNull()) {
	bool updatedSta = iEvent.getProvenance(staTrack.id()).productInstanceName() == std::string("UpdatedAtVtx");
	//h_sta_updated->Fill(0);
	if(updatedSta) h_sta_updated->Fill(1);
	else h_sta_updated->Fill(-1);
      }
      
      //reco::TrackRef staTemp;
      edm::Handle<reco::TrackToTrackMap> updatedStaAssoMap;

      if( iEvent.getByLabel(theSTACollectionLabel.label(),updatedStaAssoMap) ) {

	reco::TrackRef staUpdated;
	reco::TrackRef staNotUpdated;
	reco::TrackToTrackMap::const_iterator iEnd;
	reco::TrackToTrackMap::const_iterator iBegin;
	reco::TrackToTrackMap::const_iterator iii;

	for(iii=updatedStaAssoMap->begin(); 
	    iii!=updatedStaAssoMap->end(); ++iii) {
	  staNotUpdated = (*iii).key;
	  staUpdated = (*iii).val;

	  if(staUpdated==staTrack) {
	    LogTrace("MatchAnalyzer")<<"Checking eta flip1";
	    int etaFlip1 = 
	      ( (staUpdated->eta() * staNotUpdated->eta() ) < 0) ? -1 : 1;
	    h_sta_etaFlip1->Fill(etaFlip1);
	    LogTrace("MatchAnalyzer")<<"Checking eta flip1 " << etaFlip1;
	    if(etaFlip1 < 0) etaFliped = true;
	    //adam eta-flip-fix
	    staTrackForRegion = (fixForFlip && etaFliped) ? staNotUpdated : staUpdated;
	  }
	}
      }
    }

    //////////////////aaaaaaaaaaaaa

    if(staTrack.isAvailable() && tkTrack.isAvailable()) {
      h_muon_refit->Fill(0);
      CandidateContainer result = build(staTrack,tkTrack);
      LogTrace("MatchAnalyzer") << 
	"                        Found "<< result.size() << " GLBMuons from one STACand \n" << 
	"                           with chi2/DoF " << result.front()->trajectory()->chiSquared()/result.front()->trajectory()->ndof() << 
	"\n                           from a sta with an etaFlip? " << etaFliped;
      if(result.size()) {
	h_muon_refit->Fill(1);
	h_muon_refit_chi2->Fill(result.front()->trajectory()->chiSquared()/result.front()->trajectory()->ndof());
      }

      for( CandidateContainer::const_iterator it = result.begin(); it != result.end(); ++it) {
	if ( (*it)->trajectory() ) delete (*it)->trajectory();
	if ( (*it)->trackerTrajectory() ) delete (*it)->trackerTrajectory();
	if ( *it ) delete (*it);
      }
      result.clear();  
      

    }
    
    //////////////////////////////
    //
    // end special (temporary) selection
    //
    //////////////////////////////

    /////////////////////////////
    //
    // Make a new collection for all good / bad / all tracks
    //
    /////////////////////////////
    vector<TrackCand> tkTrackCandsNew;
    for(vector<TrackCand>::const_iterator iTk = tkTrackCands.begin();
	iTk != tkTrackCands.end(); ++iTk) {
      if(useAll==1) { 
	tkTrackCandsNew.push_back(*iTk); 
      }
      else if(useAll==2 && tkTrack.isAvailable() && iTk->second == tkTrack) {
	tkTrackCandsNew.push_back(*iTk); 
      }
      else if(useAll==3 && tkTrack.isAvailable() && iTk->second != tkTrack) {
	tkTrackCandsNew.push_back(*iTk); 
      }
    }
    
    //LogDebug("MatchAnalyzer");  
    if(iMuon->isGlobalMuon()) {
      if(glbTrack.isAvailable()) {
	glb_combined->SetPoint(iGlbMu,glbTrack->eta(),glbTrack->phi());
	glb_combined->SetPointError(iGlbMu,glbTrack->etaError(),glbTrack->phiError());
      }
      if(staTrack.isAvailable()) {
	glb_outer->SetPoint(iGlbMu,staTrack->eta(),staTrack->phi());
	glb_outer->SetPointError(iGlbMu,staTrack->etaError(),staTrack->phiError());
      }
      if(tkTrack.isAvailable()) {
	glb_inner->SetPoint(iGlbMu,tkTrack->eta(),tkTrack->phi());
	glb_inner->SetPointError(iGlbMu,tkTrack->etaError(),tkTrack->phiError());
      }
      iGlbMu++;
    }
    //LogDebug("MatchAnalyzer");  
    vector<TrackCand> tkPreCandColl;
    vector<TrackCand> tkPreCandCollFixed;

    ////////////////////////////////////
    //
    // For every muon that is a STA, check the matching observables
    //
    ////////////////////////////////////

    if(iMuon->isStandAloneMuon()) {
      TrackCand staCand((Trajectory*)(0),staTrack);     
      if(staTrack.isAvailable()) {
	sta_muon->SetPoint(iSta,staTrack->eta(),staTrack->phi());
	sta_muon->SetPointError(iSta,staTrack->etaError(),staTrack->phiError());
      }
      h_muon_steps->Fill(3);
      //LogDebug("MatchAnalyzer");  
      ///////////////////////////////////
      ///////////////////////////////////
      /*
      {
	int iTkAll = 0;
	for(vector<TrackCand>::const_iterator iTk = tkTrackCandsNew.begin();
	    iTk != tkTrackCandsNew.end(); ++iTk) {
	  std::pair<TrajectoryStateOnSurface, TrajectoryStateOnSurface> tsosPair = theTrackMatcher->convertToTSOSTk(staCand,*iTk);
	  pos_tkAllsurf->SetPoint(iTkAll,tsosPair.second.globalPosition().eta(),tsosPair.second.globalPosition().phi());
	  sta_muonsurf->SetPoint(iTkAll,tsosPair.first.globalPosition().eta(),tsosPair.first.globalPosition().phi());
	  //adam pos_tkAllsurf->SetPointError(iTkAll,tsosPair.second.globalPosition().etaError(),tsosPair.second.globalPosition().phiError());
	  iTkAll++;
	}
      }
      */
      ///////////////////////////////////
      ///////////////////////////////////

      ////////////////////////////////////
      //
      // Make collections of fixed and dynamic tracker tracks
      //
      ////////////////////////////////////
      LogTrace("MatchAnalyzer") << endl 
				<< "**********" << endl << "StaMuon " 
				<< iSta+1 << " of " << "999" << " pT " << staTrack->pt() << " eta " << staTrack->eta() << " phi " << staTrack->phi() <<endl;
      //ADAM adding the "New" for good / bad / all
      tkPreCandColl = 
	chooseRegionalTrackerTracks(staTrackForRegion,tkTrackCandsNew,iSta);
      tkPreCandCollFixed = 
	chooseRegionalTrackerTracksFixed(staTrackForRegion,tkTrackCandsNew,iSta);

      //h_nTKFixed->Fill(tkPreCandCollFixed.size());
      //h_nTKDynamic->Fill(tkPreCandColl.size());

      LogTrace("MatchAnalyzer") << "     " 
				<< endl << "Tk in Region " 
				<< tkPreCandCollFixed.size() << endl;
      h_muon_steps->Fill(4);
      /////////////////////////////
      //
      // Loop over all dynamic regional tracker tracks
      //
      /////////////////////////////
      for(vector<TrackCand>::const_iterator iTk = tkPreCandColl.begin();
	  iTk != tkPreCandColl.end(); ++iTk) {
	/////////////////////////////
	//
	// Set the point for all (dynamic) regional tracker tracks
	//
	/////////////////////////////	
	pos_tkCand->SetPoint(iTkDynamic,iTk->second->eta(),iTk->second->phi());
	pos_tkCand->SetPointError(iTkDynamic,iTk->second->etaError(),iTk->second->phiError());
	iTkDynamic++;
      }
      
      int iTkSurf = 0;
      int iiTk = 1;
      h_muon_steps->Fill(5);
      /////////////////////////////
      //
      // Loop over all fixed regional tracker tracks
      //
      /////////////////////////////
      for(vector<TrackCand>::const_iterator iTk = tkPreCandCollFixed.begin();
	  iTk != tkPreCandCollFixed.end(); ++iTk) {
	LogTrace("MatchAnalyzer") << "*****" << endl 
				  << "Tk " << iiTk << " of " 
				  << tkPreCandCollFixed.size() 
				  << " and iTkFixed " << iTkFixed 
				  << " pT " << iTk->second->pt() << endl;
	iiTk++;
	/////////////////////////////
	//
	// Set the point for all (fixed) regional tracker tracks
	//
	/////////////////////////////
	pos_tkCandFixed->SetPoint(iTkFixed,iTk->second->eta(),iTk->second->phi());
	pos_tkCandFixed->SetPointError(iTkFixed,iTk->second->etaError(),iTk->second->phiError());	
	iTkFixed++;
		
	if(etaFlipBool && iiTk==2) h_muon_steps->Fill(6);
	////////////////////////////////////
	//
	// For each surface, compute observables
	//
	////////////////////////////////////
	std::pair<TrajectoryStateOnSurface, TrajectoryStateOnSurface> tsosPair;
	std::string dirName;
	for(int i = 0; i<3; i++ ) {
	  
	  switch(i) {
	  case 0 :
	    LogTrace("MatchAnalyzer") << "     **MuHitSurface" << endl;
	    tsosPair = theTrackMatcher->convertToTSOSMuHit(staCand,*iTk);
	    dirName = std::string("matchAnalyzerMuHit");
	    LogTrace("MatchAnalyzer") << "        ConvertToMuHitSurface muon isValid " << tsosPair.first.isValid() << " tk isValid " << tsosPair.second.isValid() << endl;
	    break;
	  case 1 :
	    LogTrace("MatchAnalyzer") << "     **TkHitSurface" << endl;
	    tsosPair = theTrackMatcher->convertToTSOSTkHit(staCand,*iTk);
	    dirName = std::string("matchAnalyzerTkHit");
	    LogTrace("MatchAnalyzer") << "        ConvertToTkHitSurface muon isValid " << tsosPair.first.isValid() << " tk isValid " << tsosPair.second.isValid() << endl;
	    break;
	  case 2 :
	    LogTrace("MatchAnalyzer") << "     **TkSurface" << endl;
	    tsosPair = theTrackMatcher->convertToTSOSTk(staCand,*iTk);
	    dirName = std::string("matchAnalyzerTkSurf");
	    LogTrace("MatchAnalyzer") << "        ConvertToTkSurface muon isValid " << tsosPair.first.isValid() << " tk isValid " << tsosPair.second.isValid() << endl;
	    break;
	  default:
	    tsosPair = theTrackMatcher->convertToTSOSMuHit(staCand,*iTk);
	    LogTrace("MatchAnalyzer") << "      ConvertToMuHitSurface (default) muon isValid " << tsosPair.first.isValid() << " tk isValid " << tsosPair.second.isValid() << endl;
	  }
	  if(etaFlipBool && i==0) h_muon_steps->Fill(7);
	  std::map<std::string, TH1*> localDir = 
	    testPlots[dirName];//LogDebug("MatchAnalyzer");

	  //if(propBool) {
	  //localDir["h_propagation"]->Fill(0);
	  //if(tsosPair.first.isValid()) localDir["h_propagation"]->Fill(1);
	  //propBool = false;
	  //}
	  LogTrace("MatchAnalyzer") << "            STA before eta " 
				    << staCand.second->eta() 
				    << " phi " << staCand.second->phi() 
				    << "\n" 
				    << "                inner "<<staCand.second->innerPosition()
				    << "\n" 
				    << "                outer "<<staCand.second->outerPosition()
				    << endl;

	  LogTrace("MatchAnalyzer") << "            TK  before eta " 
				    << iTk->second->eta() << " phi " 
				    << iTk->second->phi()
				    << "\n"
				    << "                inner "<<iTk->second->innerPosition()
				    << "\n"
				    << "                outer "<<iTk->second->outerPosition()
				    << endl;
	  if(!tsosPair.first.isValid() || !tsosPair.second.isValid()) continue;
	  LogTrace("MatchAnalyzer") << "            STA after  eta " 
				    << tsosPair.first.globalPosition().eta() 
				    << " phi " << tsosPair.first.globalPosition().phi() 
				    << "\n"
				    << "                inner " << tsosPair.first.globalPosition()
				    << "\n"
				    << "                outer " << tsosPair.first.globalPosition()
				    << endl;
	  LogTrace("MatchAnalyzer") << "            TK  after  eta " 
				    << tsosPair.second.globalPosition().eta() 
				    << " phi " << tsosPair.second.globalPosition().phi() 
				    << "\n"
				    << "                inner " << tsosPair.second.globalPosition()
				    << "\n"
				    << "                outer " << tsosPair.second.globalPosition()
				    << endl;
	  
	  if(etaFlipBool && i==0) h_muon_steps->Fill(8);

	  // calculate matching variables
	  double distance = 
	    theTrackMatcher->match_d(tsosPair.first,tsosPair.second);
	  double Distance = 
	    theTrackMatcher->match_D(tsosPair.first,tsosPair.second);
	  double chi2 = 
	    theTrackMatcher->match_Chi2(tsosPair.first,tsosPair.second);
	  double loc_chi2 = 
	    theTrackMatcher->match_dist(tsosPair.first,tsosPair.second);
	  double loc_chi2_2 = 
	    theTrackMatcher->match_dist2(tsosPair.first,tsosPair.second);
	  double loc_chi2_3 = 
	    theTrackMatcher->match_dist3(tsosPair.first,tsosPair.second);
	  double deltaR = 
	    theTrackMatcher->match_Rpos(tsosPair.first,tsosPair.second);
	  
	  //LogDebug("MatchAnalyzer");
	  
	  //std::map<std::string, TH1*> localDir = testPlots[dirName];//LogDebug("MatchAnalyzer");
	  
	  localDir["h_distance"]->Fill(distance);
	  localDir["h_distanceOverStaPt"]->Fill(distance/staTrack->pt());
	  localDir["h_distanceOverStaP"]->Fill(distance/staTrack->p());
	  localDir["h_distanceOverTkPt"]->Fill(distance/iTk->second->pt());
	  localDir["h_distanceOverTkP"]->Fill(distance/iTk->second->p());
	  localDir["h_Distance"]->Fill(Distance);
	  //localDir["h_DistanceOverPt"]->Fill(Distance/staTrack->pt());
	  localDir["h_chi2"]->Fill(chi2);
	  localDir["h_loc_chi2"]->Fill(loc_chi2);
	  localDir["h_loc_chi2_2"]->Fill(loc_chi2_2);
	  localDir["h_loc_chi2_3"]->Fill(loc_chi2_3);
	  localDir["h_deltaR"]->Fill(deltaR);
	  localDir["h_distance_loc_chi2"]->Fill(distance,loc_chi2);
	  localDir["h_distance_loc_chi2_2"]->Fill(distance,loc_chi2_2);
	  localDir["h_distance_loc_chi2_3"]->Fill(distance,loc_chi2_3);	  
	  
	  
	  if(i==0) {
	    LogDebug("MatchAnalyzer") << "iTk pT " << iTk->second->pt() 
				      << " eta " << iTk->second->eta() 
				      << " phi " << iTk->second->phi() << endl; 
	    LogTrace("MatchAnalyzer") << "  distance " << distance 
				      << "    distance cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  Distance " << Distance 
				      << "    distance cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  chi2 " << chi2 
				      << "     chi2 cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  loc_chi2 " << loc_chi2 
				      << "     locChi2 cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  loc_chi2_2 " << loc_chi2_2 
				      << "     locChi2 cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  loc_chi2_3 " << loc_chi2_3 
				      << "     locChi2 cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  deltaR " << deltaR 
				      << "     deltaR cut " << " " << endl;
	    
	    LogTrace("MatchAnalyzer") << "  dR1 " << fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta()<1.5*0.2) << endl;
	    LogTrace("MatchAnalyzer") << "  dR2 " << (fabs(deltaPhi(tsosPair.second.globalPosition().phi(),tsosPair.first.globalPosition().phi())) < 0.2) << endl;
	    
	    LogTrace("MatchAnalyzer") << "  dR1redo " 
				      << fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta())
				      << (fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta() < 1.5 * 0.2) ) << " " 
				      << (fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta() ) < 1.5 * 0.2) << " " 
				      << endl;
	    
	    LogTrace("MatchAnalyzer") << "  dR2redo " << fabs(deltaPhi(tsosPair.second.globalPosition().phi(),tsosPair.first.globalPosition().phi()) ) << endl;
	 
	    double x1 = tsosPair.first.localPosition().x();
	    double y1 = tsosPair.first.localPosition().y();
	    double xx1 = tsosPair.first.localError().positionError().xx();
	    double yy1 = tsosPair.first.localError().positionError().yy();
	    
	    double x2 = tsosPair.second.localPosition().x();
	    double y2 = tsosPair.second.localPosition().y();
	    double xx2 = tsosPair.second.localError().positionError().xx();
	    double yy2 = tsosPair.second.localError().positionError().yy();
	    
	    surface_error1->SetPoint(iSta,x1,y1);
	    surface_error2->SetPoint(surfaceOffset+iTkSurf,x2,y2);
	    surface_error1->SetPointError(iSta,xx1,yy1);
	    surface_error2->SetPoint(surfaceOffset+iTkSurf,xx2,yy2);

	    
	    int etaFlip = ( (staCand.second->eta() * tsosPair.first.globalPosition().eta()) < 0) ? -1 : 1;
	    if(etaFlipBool) {
	      h_sta_etaFlip->Fill(etaFlip);
	      h_muon_steps->Fill(9);
	      etaFlipBool = false;
	    }

	    
	    iTkSurf++;
	    //tkCandCollFixed.push_back(TrackCand((Trajectory*)(0),tkRef));
	  } // end if i == 0 special prints
	  if(i==1) {
	    LogDebug("MatchAnalyzer") << "iTk pT " << iTk->second->pt() 
				      << " eta " << iTk->second->eta() 
				      << " phi " << iTk->second->phi() << endl; 
	    LogTrace("MatchAnalyzer") << "  distance " << distance 
				      << "    distance cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  Distance " << Distance 
				      << "    distance cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  chi2 " << chi2 
				      << "     chi2 cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  loc_chi2 " << loc_chi2 
				      << "     locChi2 cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  loc_chi2_2 " << loc_chi2_2 
				      << "     locChi2 cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  loc_chi2_3 " << loc_chi2_3 
				      << "     locChi2 cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  deltaR " << deltaR 
				      << "     deltaR cut " << " " << endl;
	    
	    LogTrace("MatchAnalyzer") << "  dR1 " << fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta()<1.5*0.2) << endl;
	    LogTrace("MatchAnalyzer") << "  dR2 " << (fabs(deltaPhi(tsosPair.second.globalPosition().phi(),tsosPair.first.globalPosition().phi())) < 0.2) << endl;
	    
	    LogTrace("MatchAnalyzer") << "  dR1redo " 
				      << fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta())
				      << (fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta() < 1.5 * 0.2) ) << " " 
				      << (fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta() ) < 1.5 * 0.2) << " " 
				      << endl;
	    
	    LogTrace("MatchAnalyzer") << "  dR2redo " << fabs(deltaPhi(tsosPair.second.globalPosition().phi(),tsosPair.first.globalPosition().phi()) ) << endl;
	  } // end if i == 1 special prints
	  if(i==2) {
	    LogDebug("MatchAnalyzer") << "iTk pT " << iTk->second->pt() 
				      << " eta " << iTk->second->eta() 
				      << " phi " << iTk->second->phi() << endl; 
	    LogTrace("MatchAnalyzer") << "  distance " << distance 
				      << "    distance cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  Distance " << Distance 
				      << "    distance cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  chi2 " << chi2 
				      << "     chi2 cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  loc_chi2 " << loc_chi2 
				      << "     locChi2 cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  loc_chi2_2 " << loc_chi2_2 
				      << "     locChi2 cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  loc_chi2_3 " << loc_chi2_3 
				      << "     locChi2 cut " << " " << endl;
	    LogTrace("MatchAnalyzer") << "  deltaR " << deltaR 
				      << "     deltaR cut " << " " << endl;
	    
	    LogTrace("MatchAnalyzer") << "  dR1 " << fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta()<1.5*0.2) << endl;
	    LogTrace("MatchAnalyzer") << "  dR2 " << (fabs(deltaPhi(tsosPair.second.globalPosition().phi(),tsosPair.first.globalPosition().phi())) < 0.2) << endl;
	    
	    LogTrace("MatchAnalyzer") << "  dR1redo " 
				      << fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta())
				      << (fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta() < 1.5 * 0.2) ) << " " 
				      << (fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta() ) < 1.5 * 0.2) << " " 
				      << endl;
	    
	    LogTrace("MatchAnalyzer") << "  dR2redo " << fabs(deltaPhi(tsosPair.second.globalPosition().phi(),tsosPair.first.globalPosition().phi()) ) << endl;
	 
	  } // end if i == 2 special prints	  
	} // end loop over each surface	
      } // end loop over all fixed regional tracker tracks

      surfaceOffset = surfaceOffset + iTkSurf + 1;

      ////////////////////////////////////////
      //
      // Run the general matcher sta to track collection
      //
      ////////////////////////////////////////
      
      /* ADAM: let's call the general matcher only once per analyzer
	 vector<TrackCand> selectedTrackerTracks = theTrackMatcher->match(TrackCand((Trajectory*)(0),staTrack), tkPreCandColl);
	 for(vector<TrackCand>::const_iterator iTk=selectedTrackerTracks.begin();
	 iTk != selectedTrackerTracks.end(); ++iTk) {
	 pos_selectedTkCand->SetPoint(iSelTkDynamic,iTk->second->eta(),iTk->second->phi());
	 pos_selectedTkCand->SetPointError(iSelTkDynamic,iTk->second->etaError(),iTk->second->phiError());
	 iSelTkDynamic++;
	 }
      */
      
      vector<TrackCand> selectedTrackerTracksFixed = theTrackMatcher->match(TrackCand((Trajectory*)(0),staTrack), tkPreCandCollFixed);
      //adam: try the local distance crap here
      //end adam
      //successfull match here
      
      if(selectedTrackerTracksFixed.size()) h_muon_steps->Fill(10);

      for(vector<TrackCand>::const_iterator iTk=selectedTrackerTracksFixed.begin();
	  iTk != selectedTrackerTracksFixed.end(); ++iTk) {
	LogTrace("MatchAnalyzer") << "-----" << endl 
				  << "selected pt " << iTk->second->pt() 
				  << " eta " << iTk->second->eta() 
				  << " phi " << iTk->second->phi() << endl; 
	pos_selectedTkCandFixed->SetPoint(iSelTkFixed,iTk->second->eta(),iTk->second->phi());
	pos_selectedTkCandFixed->SetPointError(iSelTkFixed,iTk->second->etaError(),iTk->second->phiError());
	iSelTkFixed++;
      }
      
      iSta++;
    } //end loop over STA

    ////////////////////
    //
    // Set the Points for all tracker muons
    //
    ////////////////////

    if(iMuon->isTrackerMuon()) {
      if(tkTrack.isAvailable()) {
	tk_muon->SetPoint(iTkMu,tkTrack->eta(),tkTrack->phi());
	tk_muon->SetPointError(iTkMu,tkTrack->etaError(),tkTrack->phiError());
      }
      iTkMu++;
    }

  } // end loop over all muons
}


// ------------ method called once each job just before starting event loop  ------------
void 
GlobalMatchingAnalyser::beginJob()
{
  edm::Service<TFileService> fs;


  TFileDirectory subDir = fs->mkdir( "matchAnalyzer" );

  ////tg1 = new TGraphErrors();
  //tg1 = subDir.make<TGraphErrors>(10); // TGraphErrors();
  //tg1->SetName("tg1_name");
  //tg1->SetTitle("tg1_title");
  //tg1->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  //tg1->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  glb_combined = subDir.make<TGraphErrors>();
  glb_combined->SetName("glb_combined");
  glb_combined->SetTitle("Global Combined Muon");
  glb_combined->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  glb_combined->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  glb_combined->SetLineColor(4);

  glb_inner = subDir.make<TGraphErrors>();
  glb_inner->SetName("glb_inner");
  glb_inner->SetTitle("Global Inner Muon");
  glb_inner->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  glb_inner->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  glb_inner->SetLineColor(3);

  glb_outer = subDir.make<TGraphErrors>();
  glb_outer->SetName("glb_outer");
  glb_outer->SetTitle("Global Outer Muon");
  glb_outer->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  glb_outer->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  glb_outer->SetLineColor(2);

  sta_muon = subDir.make<TGraphErrors>();
  sta_muon->SetName("sta_muon");
  sta_muon->SetTitle("Stand-alone Muon");
  sta_muon->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  sta_muon->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  sta_muon->SetLineColor(kRed);
  sta_muon->SetMarkerStyle(3);
  sta_muon->SetMarkerColor(2);

  sta_muonsurf = subDir.make<TGraphErrors>();
  sta_muonsurf->SetName("sta_muonsurf");
  sta_muonsurf->SetTitle("Stand-alone Muon");
  sta_muonsurf->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  sta_muonsurf->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  sta_muonsurf->SetLineColor(kOrange);
  sta_muonsurf->SetMarkerStyle(3);
  sta_muonsurf->SetMarkerColor(2);

  tk_muon = subDir.make<TGraphErrors>();
  tk_muon->SetName("tk_muon");
  tk_muon->SetTitle("Tracker Muon");
  tk_muon->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  tk_muon->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  tk_muon->SetMarkerStyle(28);
  tk_muon->SetMarkerColor(3);

  region_fixed = subDir.make<TGraphErrors>();
  region_fixed->SetName("region_fixed");
  region_fixed->SetTitle("Fixed Region Size");
  region_fixed->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  region_fixed->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  region_dynamic = subDir.make<TGraphErrors>();
  region_dynamic->SetName("region_dynamic");
  region_dynamic->SetTitle("Dynamic Region Size");
  region_dynamic->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  region_dynamic->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  pos_tkAll = subDir.make<TGraphErrors>();
  pos_tkAll->SetName("pos_tkAll");
  pos_tkAll->SetTitle("All TK Tracks");
  pos_tkAll->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_tkAll->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  pos_tkAll->GetHistogram()->GetXaxis()->SetTitle("#eta");
  pos_tkAll->GetHistogram()->GetYaxis()->SetTitle("#phi");
  pos_tkAll->SetLineColor(kBlue);

  pos_tkAllsurf = subDir.make<TGraphErrors>();
  pos_tkAllsurf->SetName("pos_tkAllsurf");
  pos_tkAllsurf->SetTitle("All TK Tracks on Tracker Surface");
  pos_tkAllsurf->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_tkAllsurf->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  pos_tkAllsurf->GetHistogram()->GetXaxis()->SetTitle("#eta");
  pos_tkAllsurf->GetHistogram()->GetYaxis()->SetTitle("#phi");
  pos_tkAllsurf->SetLineColor(kOrange);

  pos_tkCand = subDir.make<TGraphErrors>();
  pos_tkCand->SetName("pos_tkCand");
  pos_tkCand->SetTitle("TK Candidates");
  pos_tkCand->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_tkCand->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  pos_tkCandFixed = subDir.make<TGraphErrors>();
  pos_tkCandFixed->SetName("pos_tkCandFixed");
  pos_tkCandFixed->SetTitle("TK Candidates");
  pos_tkCandFixed->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_tkCandFixed->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  pos_tkCandFixed->SetMarkerStyle(25);

  pos_selectedTkCand = subDir.make<TGraphErrors>();
  pos_selectedTkCand->SetName("pos_selectedTkCand");
  pos_selectedTkCand->SetTitle("Matched TK Candidates");
  pos_selectedTkCand->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_selectedTkCand->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  pos_selectedTkCand->SetMarkerStyle(24);

  pos_selectedTkCandFixed = subDir.make<TGraphErrors>();
  pos_selectedTkCandFixed->SetName("pos_selectedTkCandFixed");
  pos_selectedTkCandFixed->SetTitle("Matched TK Candidates");
  pos_selectedTkCandFixed->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_selectedTkCandFixed->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  pos_selectedTkCandFixed->SetMarkerStyle(24);
  
  surface_error1 = subDir.make<TGraphErrors>();
  surface_error1->SetName("surface_error1");
  surface_error1->SetTitle("Common Surface1");
  //surface_error1->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  //surface_error1->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  surface_error2 = subDir.make<TGraphErrors>();
  surface_error2->SetName("surface_error2");
  surface_error2->SetTitle("Common Surface2");
  //surface_error2->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  //surface_error2->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  //TFileDirectory testSubDir = fs->mkdir( "testMatchAnalyzer" );
  //(testPlots["testMatchAnalyzer"])["h_distance_test"] = testSubDir.make<TH1F>("h_distance_test","distance",100,0,50);

  h_sta_pt = subDir.make<TH1F>("h_sta_pt","STA p_{T}",100,0.,50.);
  h_sta_p = subDir.make<TH1F>("h_sta_p","STA p",100,0.,50.);
  h_sta_rho = subDir.make<TH1F>("h_sta_rho","STA #rho",100,0.,50.);
  h_sta_R = subDir.make<TH1F>("h_sta_R","STA R",100,0.,50.);
  h_sta_rho_outer = subDir.make<TH1F>("h_sta_rho_outer","STA #rho_{outer}",100,0.,50.);
  h_sta_R_outer = subDir.make<TH1F>("h_sta_R_outer","STA R_{outer}",100,0.,50.);
  h_sta_cut = subDir.make<TH1F>("h_sta_cut","STA Pass Cut",5,-2.5,2.5);
  h_sta_etaFlip = subDir.make<TH1F>("h_sta_etaFlip","#eta_{STAbefore} * #eta_{STAfter}",3,-1.5,1.5);
  h_sta_etaFlip1 = subDir.make<TH1F>("h_sta_etaFlip1","#eta_{STAupdated} * #eta_{STAnotUpdated}",3,-1.5,1.5);
  h_sta_updated = subDir.make<TH1F>("h_sta_updated","STA was UpdatedAtVtx",3,-1.5,1.5);
  h_muon_steps = subDir.make<TH1F>("h_muon_steps","Muon Steps",11,-0.5,10.5);
  h_muon_refit = subDir.make<TH1F>("h_muon_refit","Refit Muon",3,-1.5,1.5);
  h_muon_refit_chi2 = subDir.make<TH1F>("h_muon_refit_chi2","Refit Muon #chi^2",100,0.,100.);
  h_muon_type = subDir.make<TH1F>("h_muon_type","Muon isGLB, isTK, isSTA",4,-0.5,3.5);

  std::vector<std::string> dirName;
  dirName.push_back("matchAnalyzerMuHit");
  dirName.push_back("matchAnalyzerTkHit");
  dirName.push_back("matchAnalyzerTkSurf");
  
  for(std::vector<std::string>::const_iterator iDir=dirName.begin(); iDir != dirName.end(); ++iDir) {
    
    TFileDirectory subDir2 = fs->mkdir( *iDir );
    
    //(testPlots[dirName.front()])["h_distance_muHitNew"] = subDir2.make<TH1F>("h_distancetest","distancetest",100,0,50);
    
    (testPlots[*iDir])["h_distance"] = subDir2.make<TH1F>("h_distance","distance",100,0,50);
    (testPlots[*iDir])["h_distanceOverStaPt"] = subDir2.make<TH1F>("h_distanceOverStaPt","distance/p_{T}^{STA}",100,0,5);
    (testPlots[*iDir])["h_distanceOverStaP"] = subDir2.make<TH1F>("h_distanceOverStaP","distance/p^{STA}",100,0,5);
    (testPlots[*iDir])["h_distanceOverTkPt"] = subDir2.make<TH1F>("h_distanceOverTkPt","distance/p_{T}^{TK}",100,0,5);
    (testPlots[*iDir])["h_distanceOverTkP"] = subDir2.make<TH1F>("h_distanceOverTkP","distance/p^{TK}",100,0,5);
    (testPlots[*iDir])["h_Distance"] = subDir2.make<TH1F>("h_Distance","Distance",100,0,50);
    //(testPlots[*iDir])["h_DistanceOverPt"] = subDir2.make<TH1F>("h_DistanceOverPt","Distance/p_{T}",100,0,5);
    (testPlots[*iDir])["h_chi2"] = subDir2.make<TH1F>("h_chi2","chi2",100,0,500);
    (testPlots[*iDir])["h_loc_chi2"] = subDir2.make<TH1F>("h_loc_chi2","loc_chi2",100,0,0.001);
    (testPlots[*iDir])["h_loc_chi2_2"] = subDir2.make<TH1F>("h_loc_chi2_2","loc_chi2_2",500,0,500.);
    (testPlots[*iDir])["h_loc_chi2_3"] = subDir2.make<TH1F>("h_loc_chi2_3","loc_chi2_3",500,0,500.);
    (testPlots[*iDir])["h_deltaR"] = subDir2.make<TH1F>("h_deltaR","deltaR",100,0,2.);
    (testPlots[*iDir])["h_distance_loc_chi2"] = subDir2.make<TH2F>("h_distance_loc_chi2"," loc_chi2 vs distance",100,0,50,100,0,500.);
    (testPlots[*iDir])["h_distance_loc_chi2_2"] = subDir2.make<TH2F>("h_distance_loc_chi2_2"," loc_chi2_2 vs distance",100,0,50,100,0,500.);
    (testPlots[*iDir])["h_distance_loc_chi2_3"] = subDir2.make<TH2F>("h_distance_loc_chi2_3"," loc_chi2_3 vs distance",100,0,50,100,0,500.);

    (testPlots[*iDir])["h_propagation"] = subDir2.make<TH1F>("h_propagation","Propagation to same surface",3,-1.5,1.5);
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
GlobalMatchingAnalyser::endJob() {
  
  
}

//
// define a region of interest within the tracker
//
RectangularEtaPhiTrackingRegion 
GlobalMatchingAnalyser::defineRegionOfInterest(const reco::TrackRef& staTrack) const {

  RectangularEtaPhiTrackingRegion* region1 = theRegionBuilder->region(staTrack);
  
  TkTrackingRegionsMargin<float> etaMargin(fabs(region1->etaRange().min() - region1->etaRange().mean()),
					   fabs(region1->etaRange().max() - region1->etaRange().mean()));
  
  RectangularEtaPhiTrackingRegion region2(region1->direction(),
					  region1->origin(),
					  region1->ptMin(),
					  region1->originRBound(),
					  region1->originZBound(),
					  etaMargin,
					  region1->phiMargin());
  
  delete region1;
  return region2;
  
}

//
// select tracker tracks within a region of interest
//
vector<TrackCand> 
GlobalMatchingAnalyser::chooseRegionalTrackerTracksFixed(const reco::TrackRef& staCand, 
                                                         const vector<TrackCand>& tkTs,
							 int iSta) {
  
  // define eta-phi region
  RectangularEtaPhiTrackingRegion regionOfInterest = defineRegionOfInterest(staCand);
  
  // get region's etaRange and phiMargin
  PixelRecoRange<float> etaRange = regionOfInterest.etaRange();
  TkTrackingRegionsMargin<float> phiMargin = regionOfInterest.phiMargin();

  region_fixed->SetPoint(iSta,regionOfInterest.direction().eta(),regionOfInterest.direction().phi());
  region_fixed->SetPointError(iSta,1.,1.);

  vector<TrackCand> result;

  double deltaR_max = 1.0;

  for ( vector<TrackCand>::const_iterator is = tkTs.begin(); 
	is != tkTs.end(); ++is ) {

    double deltaR_tmp = deltaR(static_cast<double>(regionOfInterest.direction().eta()),
			       static_cast<double>(regionOfInterest.direction().phi()),
			       is->second->eta(), is->second->phi());

    // for each trackCand in region, add trajectory and add to result
    if (deltaR_tmp < deltaR_max) {
      TrackCand tmpCand = TrackCand(*is);
      result.push_back(tmpCand);
    }
  }

  return result; 

}

//
// select tracker tracks within a region of interest
//
vector<TrackCand> 
GlobalMatchingAnalyser::chooseRegionalTrackerTracks( const TrackRef& staCand, 
						     const vector<TrackCand>& tkTs,
						     int iSta) {
  
  // define eta-phi region
  RectangularEtaPhiTrackingRegion region = defineRegionOfInterest(staCand);
  
  //if the limits are too stringent rescale limits
  typedef PixelRecoRange< float > Range;
  typedef TkTrackingRegionsMargin< float > Margin;
  
  Range etaRange   = region.etaRange();
  double etaLimit  = (fabs(fabs(etaRange.max()) - fabs(etaRange.mean())) <0.1) ? 0.1 : fabs(fabs(etaRange.max()) - fabs(etaRange.mean())) ;
  
  Margin phiMargin = region.phiMargin();
  double phiLimit  = (phiMargin.right() < 0.1 ) ? 0.1 : phiMargin.right(); 
  
  region_dynamic->SetPoint(iSta,region.direction().eta(),region.direction().phi());
  region_dynamic->SetPointError(iSta,etaLimit,phiLimit);

  vector<TrackCand> result;

  for ( vector<TrackCand>::const_iterator is = tkTs.begin(); 
	is != tkTs.end(); ++is ) {
    
    double trackEta = is->second->eta();
    double trackPhi = is->second->phi();
    
    // Clean  
    bool inEtaRange = trackEta >= (etaRange.mean() - etaLimit) && trackEta <= (etaRange.mean() + etaLimit) ;
    bool inPhiRange = (fabs(deltaPhi(trackPhi,double(region.direction().phi()))) < phiLimit );
    
    if(inEtaRange && inPhiRange) result.push_back(*is);
    
  }

  return result;

}


GlobalMatchingAnalyser::CandidateContainer GlobalMatchingAnalyser::build(const TrackRef staCand, const TrackRef tkCand) const {

  LogTrace("MatchAnalyzer")<<
    "                     ***Begin Build";

  CandidateContainer refittedResult;
  ConstRecHitContainer muonRecHits = getTransientRecHits(*(staCand));

  // check order of muon measurements
  if ( (muonRecHits.size() > 1) &&
       ( muonRecHits.front()->globalPosition().mag() >
	 muonRecHits.back()->globalPosition().mag() ) ) {
    LogTrace("MuonAnalyzer")<< "   reverse order: ";
    stable_sort(muonRecHits.begin(),muonRecHits.end(),RecHitLessByDet(alongMomentum));
  }
  
  stable_sort(muonRecHits.begin(),muonRecHits.end(),ComparatorInOut());
  
  // cut on tracks with low momenta
  LogTrace("MuonAnalyzer")<< "   Track p and pT " << tkCand->p() << " " << tkCand->pt();
  //if(  tkCand->p() < thePCut || tkCand->pt() < thePtCut  ) continue;
  
  ConstRecHitContainer trackerRecHits;
  if (tkCand.isNonnull()) {
    trackerRecHits = getTransientRecHits(*tkCand);
  } else {
    LogDebug("MuonAnalyzer")<<"     NEED HITS FROM TRAJ";
    //trackerRecHits = (*it)->trackerTrajectory()->recHits();
  }
  
  // check for single TEC RecHits in trajectories in the overalp region
  /*
    if ( fabs(tkCand->eta()) > 0.95 && fabs(tkCand->eta()) < 1.15 && tkCand->pt() < 60 ) {
    if ( theTECxScale < 0 || theTECyScale < 0 )
    trackerRecHits = selectTrackerHits(trackerRecHits);
    else
    fixTEC(trackerRecHits,theTECxScale,theTECyScale);
    }
  */
  
  RefitDirection recHitDir = checkRecHitsOrdering(trackerRecHits);
  if ( recHitDir == outToIn ) reverse(trackerRecHits.begin(),trackerRecHits.end());
  
  reco::TransientTrack tTT(tkCand,&*theService->magneticField(),theService->trackingGeometry());
  TrajectoryStateOnSurface innerTsos = tTT.innermostMeasurementState();
  
  edm::RefToBase<TrajectorySeed> tmpSeed;
  if(tkCand->seedRef().isAvailable()) tmpSeed = tkCand->seedRef();
  
  if ( !innerTsos.isValid() ) {
    LogTrace("MuonAnalyzer") << "     inner Trajectory State is invalid. ";
    return refittedResult;
  }
  
  innerTsos.rescaleError(100.);
  
  TC refitted0,refitted1;
  MuonCandidate* finalTrajectory = 0;
  Trajectory *tkTrajectory = 0;
  
  // tracker only track
  if ( true ) /*|| ! ((*it)->trackerTrajectory() 
		&& (*it)->trackerTrajectory()->isValid()) ) */ { 
    refitted0 = theTrackTransformer->transform(tkCand) ;
    if (!refitted0.empty()) tkTrajectory 
      = new Trajectory(*(refitted0.begin())); 
    else LogWarning("MatchAnalyser")
      << "     Failed to load tracker track trajectory";
  } else tkTrajectory = 0;
  if (tkTrajectory) tkTrajectory->setSeedRef(tmpSeed);
  
  // full track with all muon hits using theGlbRefitter    
  ConstRecHitContainer allRecHits = trackerRecHits;
  allRecHits.insert(allRecHits.end(), muonRecHits.begin(),muonRecHits.end());
  refitted1 = theGlbRefitter->refit( *tkCand, tTT, allRecHits,theMuonHitsOption);
  LogTrace("MatchAnalyser")<<"     This track-sta refitted to " << refitted1.size() << " trajectories";
  
  Trajectory *glbTrajectory1 = 0;
  if (!refitted1.empty()) glbTrajectory1 = new Trajectory(*(refitted1.begin()));
  else LogDebug("MatchAnalyser")<< "     Failed to load global track trajectory 1"; 
  if (glbTrajectory1) glbTrajectory1->setSeedRef(tmpSeed);
  
  finalTrajectory = 0;
  if(glbTrajectory1 && tkTrajectory) finalTrajectory 
    = new MuonCandidate(glbTrajectory1, staCand, tkCand, 
			tkTrajectory? new Trajectory(*tkTrajectory) : 0);
  
  if ( finalTrajectory ) 
    refittedResult.push_back(finalTrajectory);
  
  if(tkTrajectory) delete tkTrajectory;
  
  // choose the best global fit for this Standalone Muon based on the track probability
  CandidateContainer selectedResult;
  MuonCandidate* tmpCand = 0;
  if ( refittedResult.size() > 0 ) tmpCand = *(refittedResult.begin());
  double minProb = 9999;

  for (CandidateContainer::const_iterator iter=refittedResult.begin(); iter != refittedResult.end(); iter++) {
    double prob = trackProbability(*(*iter)->trajectory());
    LogTrace("MatchAnalyser")<<"   refitted-track-sta with pT " << (*iter)->trackerTrack()->pt() << " has probability " << prob;

    if (prob < minProb) {
      minProb = prob;
      tmpCand = (*iter);
    }
  }

  if ( tmpCand )  selectedResult.push_back(new MuonCandidate(new Trajectory(*(tmpCand->trajectory())), tmpCand->muonTrack(), tmpCand->trackerTrack(), 
							     (tmpCand->trackerTrajectory())? new Trajectory( *(tmpCand->trackerTrajectory()) ):0 ) );

  for (CandidateContainer::const_iterator it = refittedResult.begin(); it != refittedResult.end(); ++it) {
    if ( (*it)->trajectory() ) delete (*it)->trajectory();
    if ( (*it)->trackerTrajectory() ) delete (*it)->trackerTrajectory();
    if ( *it ) delete (*it);
  }
  refittedResult.clear();

  return selectedResult;


}

//
// get transient RecHits
//
TransientTrackingRecHit::ConstRecHitContainer
GlobalMatchingAnalyser::getTransientRecHits(const reco::Track& track) const {

  TransientTrackingRecHit::ConstRecHitContainer result;
  
  TrajectoryStateTransform tsTrans;

  TrajectoryStateOnSurface currTsos = tsTrans.innerStateOnSurface(track, *theService->trackingGeometry(), &*theService->magneticField());

  for (trackingRecHit_iterator hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
    if((*hit)->isValid()) {
      DetId recoid = (*hit)->geographicalId();
      if ( recoid.det() == DetId::Tracker ) {
	TransientTrackingRecHit::RecHitPointer ttrhit = theTrackerRecHitBuilder->build(&**hit);
	TrajectoryStateOnSurface predTsos =  theService->propagator(theTrackerPropagatorName)->propagate(currTsos, theService->trackingGeometry()->idToDet(recoid)->surface());

	if ( !predTsos.isValid() ) {
	  edm::LogError("MissingTransientHit")
	    <<"Could not get a tsos on the hit surface. We will miss a tracking hit.";
	  continue; 
	}
	currTsos = predTsos;
	TransientTrackingRecHit::RecHitPointer preciseHit = ttrhit->clone(predTsos);
	result.push_back(preciseHit);
      } else if ( recoid.det() == DetId::Muon ) {
	if ( (*hit)->geographicalId().subdetId() == 3 && !theRPCInTheFit) {
	  LogDebug("MatchAnalyzer") << "RPC Rec Hit discarded"; 
	  continue;
	}
	result.push_back(theMuonRecHitBuilder->build(&**hit));
      }
    }
  }
  
  return result;
}

//
// check order of RechIts on a trajectory
//
GlobalMatchingAnalyser::RefitDirection
GlobalMatchingAnalyser::checkRecHitsOrdering(const TransientTrackingRecHit::ConstRecHitContainer& recHits) const {

  if ( !recHits.empty() ) {
    ConstRecHitContainer::const_iterator frontHit = recHits.begin();
    ConstRecHitContainer::const_iterator backHit  = recHits.end() - 1;
    while ( !(*frontHit)->isValid() && frontHit != backHit ) {frontHit++;}
    while ( !(*backHit)->isValid() && backHit != frontHit )  {backHit--;}

    double rFirst = (*frontHit)->globalPosition().mag();
    double rLast  = (*backHit) ->globalPosition().mag();

    if ( rFirst < rLast ) return inToOut;
    else if (rFirst > rLast) return outToIn;
    else {
      LogError("MatchAnalyzer") << "Impossible to determine the rechits order" << endl;
      return undetermined;
    }
  }
  else {
    LogError("MatchAnalyzer") << "Impossible to determine the rechits order" << endl;
    return undetermined;
  }
}

//
// calculate the tail probability (-ln(P)) of a fit
//
double 
GlobalMatchingAnalyser::trackProbability(const Trajectory& track) const {

  if ( track.ndof() > 0 && track.chiSquared() > 0 ) { 
    return -LnChiSquaredProbability(track.chiSquared(), track.ndof());
  } else {
    return 0.0;
  }

}


//define this as a plug-in
DEFINE_FWK_MODULE(GlobalMatchingAnalyser);
