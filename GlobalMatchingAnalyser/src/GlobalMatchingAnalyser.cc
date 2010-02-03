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
// $Id: GlobalMatchingAnalyser.cc,v 1.8 2010/02/02 21:56:15 aeverett Exp $
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

#include "TFile.h"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>

//
// class decleration
//

class GlobalMatchingAnalyser : public edm::EDAnalyzer {
public:
  explicit GlobalMatchingAnalyser(const edm::ParameterSet&);
  ~GlobalMatchingAnalyser();

  typedef std::pair<const Trajectory*, reco::TrackRef> TrackCand;
  
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
  
  // ----------member data ---------------------------
  MuonTrackingRegionBuilder* theRegionBuilder;

  MuonServiceProxy* theService;

  GlobalMuonTrackMatcher* theTrackMatcher;

  TGraphErrors *tg1;
  TGraphErrors *glb_combined, *glb_inner, *glb_outer;
  TGraphErrors *sta_muon;
  TGraphErrors *tk_muon;

  TGraphErrors *region_dynamic, *region_fixed;

  TGraphErrors *pos_tkCand, *pos_tkCandFixed, *pos_tkAll;
  TGraphErrors *pos_selectedTkCand, *pos_selectedTkCandFixed;

  TGraphErrors *surface_error1;
  TGraphErrors *surface_error2;

  edm::InputTag theTrackLabel;
  edm::InputTag theMuonLabel;

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

  // the services
  ParameterSet serviceParameters = iConfig.getParameter<ParameterSet>("ServiceParameters");  
  theService = new MuonServiceProxy(serviceParameters);
  
  // the region builder
  ParameterSet regionBuilderPSet = iConfig.getParameter<ParameterSet>("MuonTrackingRegionBuilder");  
  theRegionBuilder = new MuonTrackingRegionBuilder(regionBuilderPSet,theService);

  ParameterSet trackMatcherPSet = iConfig.getParameter<ParameterSet>("GlobalMuonTrackMatcher");
  theTrackMatcher = new GlobalMuonTrackMatcher(trackMatcherPSet,theService);

}


GlobalMatchingAnalyser::~GlobalMatchingAnalyser()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  if(theRegionBuilder) delete theRegionBuilder;
  if(theService) delete theService;
  if(theTrackMatcher) delete theTrackMatcher;

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


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

   LogTrace("MatchAnalyzer") << "********************" << "Run " << iEvent.id().run() << " Event " << iEvent.id().event() ;

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

  vector<TrackCand> tkTrackCands;
  
  for ( unsigned int position = 0; position != allTrackerTracks->size(); ++position ) {
    reco::TrackRef tkTrackRef(allTrackerTracks,position);
    TrackCand tkCand = TrackCand((Trajectory*)(0),tkTrackRef);
    tkTrackCands.push_back(tkCand); 
  }

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
  
  tg1->SetPoint(0,1,1);
  tg1->SetPoint(1,2,1);
  tg1->SetPointError(0,1,.5);
  tg1->SetPointError(1,0.25,1.);

  int nMu = muonColl.size();

  int iMu = 0;
  int iSta = 0;
  int surfaceOffset = 0;
  int iTkFixed = 0;
  int iTkDynamic = 0;
  int iSelTkFixed = 0;
  int iSelTkDynamic = 0;
  for(View<Muon>::const_iterator iMuon = muonColl.begin();
      iMuon != muonColl.end(); ++iMuon, iMu++) {
    LogTrace("MatchAnalyzer") << "*****" << endl << "Muon " << iMu+1 << " of " << nMu << endl;
    const reco::TrackRef glbTrack = ( iMuon->isGlobalMuon()) ? 
      iMuon->combinedMuon() : reco::TrackRef();

    const reco::TrackRef tkTrack = ( iMuon->isTrackerMuon() ) ? 
      iMuon->innerTrack() : TrackRef();
    
    const reco::TrackRef staTrack = ( iMuon->isStandAloneMuon() ) ? 
      iMuon->outerTrack() : TrackRef();

    if(iMuon->isGlobalMuon()) {
      if(glbTrack.isAvailable()) {
	glb_combined->SetPoint(iMu,glbTrack->eta(),glbTrack->phi());
	glb_combined->SetPointError(iMu,glbTrack->etaError(),glbTrack->phiError());
      }
      if(staTrack.isAvailable()) {
	glb_outer->SetPoint(iMu,staTrack->eta(),staTrack->phi());
	glb_outer->SetPointError(iMu,staTrack->etaError(),staTrack->phiError());
      }
      if(tkTrack.isAvailable()) {
	glb_inner->SetPoint(iMu,tkTrack->eta(),tkTrack->phi());
	glb_inner->SetPointError(iMu,tkTrack->etaError(),tkTrack->phiError());
      }
    }

    vector<TrackCand> tkPreCandColl;
    vector<TrackCand> tkPreCandCollFixed;

    if(iMuon->isStandAloneMuon()) {
      TrackCand staCand((Trajectory*)(0),staTrack);     
      if(staTrack.isAvailable()) {
	sta_muon->SetPoint(iMu,staTrack->eta(),staTrack->phi());
	sta_muon->SetPointError(iMu,staTrack->etaError(),staTrack->phiError());
      }

      LogTrace("MatchAnalyzer") << endl << "**********" << endl << "StaMuon " << iSta+1 << " of " << "999" << endl;
      
      tkPreCandColl = chooseRegionalTrackerTracks(staTrack,tkTrackCands,iSta);
      tkPreCandCollFixed = chooseRegionalTrackerTracksFixed(staTrack,tkTrackCands,iSta);

      //h_nTKFixed->Fill(tkPreCandCollFixed.size());
      //h_nTKDynamic->Fill(tkPreCandColl.size());

      LogTrace("MatchAnalyzer") << "     " << endl << "Tk in Region " << tkPreCandCollFixed.size() << endl;



      //vector<TrackCand> tkCandColl;
      //vector<TrackCand> tkCandCollFixed;
      //int position = 0;
      for(vector<TrackCand>::const_iterator iTk = tkPreCandColl.begin();
	  iTk != tkPreCandColl.end(); ++iTk) {
	//TrackRef tkRef(tkPreCandColl,position);
	//position++;
	pos_tkCand->SetPoint(iTkDynamic,iTk->second->eta(),iTk->second->phi());
	pos_tkCand->SetPointError(iTkDynamic,iTk->second->etaError(),iTk->second->phiError());
	iTkDynamic++;

	//tkCandColl.push_back(TrackCand((Trajectory*)(0),tkRef));
      }

      //position = 0;
      int iTkSurf = 0;
      int iiTk = 1;
      for(vector<TrackCand>::const_iterator iTk = tkPreCandCollFixed.begin();
	  iTk != tkPreCandCollFixed.end(); ++iTk) {
	LogTrace("MatchAnalyzer") << "*****" << endl << "Tk " << iiTk << " of " << tkPreCandCollFixed.size() << endl;
	iiTk++;
	//TrackRef tkRef(tkPreCandCollFixed,position);
	//position++;
	pos_tkCandFixed->SetPoint(iTkFixed,iTk->second->eta(),iTk->second->phi());
	pos_tkCandFixed->SetPointError(iTkFixed,iTk->second->etaError(),iTk->second->phiError());

	std::pair<TrajectoryStateOnSurface, TrajectoryStateOnSurface>
	  tsosPair = theTrackMatcher->convertToTSOSMuHit(staCand,*iTk);

	LogTrace("MatchAnalyzer") << "ConvertToMuHitSurface muon isValid " << tsosPair.first.isValid() << " tk isValid " << tsosPair.second.isValid() << endl;

	if(!tsosPair.first.isValid() || !tsosPair.second.isValid()) continue;

	// calculate matching variables
	double distance = theTrackMatcher->match_d(tsosPair.first,tsosPair.second);
	double chi2 = theTrackMatcher->match_Chi2(tsosPair.first,tsosPair.second);
	double loc_chi2 = theTrackMatcher->match_dist(tsosPair.first,tsosPair.second);
	double deltaR = theTrackMatcher->match_Rpos(tsosPair.first,tsosPair.second);
	LogTrace("MatchAnalyzer") << "iTk pT " << iTk->second->pt() << " eta " << iTk->second->eta() << " phi " << iTk->second->phi() << endl; 
	LogTrace("MatchAnalyzer") << "  distance " << distance << "    distance cut " << " " << endl;
	LogTrace("MatchAnalyzer") << "  chi2 " << chi2 << "     chi2 cut " << " " << endl;
	LogTrace("MatchAnalyzer") << "  loc_chi2 " << loc_chi2 << "     locChi2 cut " << " " << endl;
	LogTrace("MatchAnalyzer") << "  deltaR " << deltaR << "     deltaR cut " << " " << endl;

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
	
	iTkSurf++;
	iTkFixed++;

	//tkCandCollFixed.push_back(TrackCand((Trajectory*)(0),tkRef));
      }

      surfaceOffset = surfaceOffset + iTkSurf + 1;

      vector<TrackCand> selectedTrackerTracks = theTrackMatcher->match(TrackCand((Trajectory*)(0),staTrack), tkPreCandColl);
      for(vector<TrackCand>::const_iterator iTk=selectedTrackerTracks.begin();
	  iTk != selectedTrackerTracks.end(); ++iTk) {
	pos_selectedTkCand->SetPoint(iSelTkDynamic,iTk->second->eta(),iTk->second->phi());
	pos_selectedTkCand->SetPointError(iSelTkDynamic,iTk->second->etaError(),iTk->second->phiError());
	iSelTkDynamic++;
      }


      vector<TrackCand> selectedTrackerTracksFixed = theTrackMatcher->match(TrackCand((Trajectory*)(0),staTrack), tkPreCandCollFixed);
      for(vector<TrackCand>::const_iterator iTk=selectedTrackerTracksFixed.begin();
	  iTk != selectedTrackerTracksFixed.end(); ++iTk) {
	LogTrace("MatchAnalyzer") << "-----" << endl << "selected pt " << iTk->second->pt() << " eta " << iTk->second->eta() << " phi " << iTk->second->phi() << endl; 
	pos_selectedTkCandFixed->SetPoint(iSelTkFixed,iTk->second->eta(),iTk->second->phi());
	pos_selectedTkCandFixed->SetPointError(iSelTkFixed,iTk->second->etaError(),iTk->second->phiError());
	iSelTkFixed++;
      }
      
      iSta++;
    } //end loop over STA

    if(iMuon->isTrackerMuon()) {
      if(tkTrack.isAvailable()) {
	tk_muon->SetPoint(iMu,tkTrack->eta(),tkTrack->phi());
	tk_muon->SetPointError(iMu,tkTrack->etaError(),tkTrack->phiError());
      }
    }


    //iMu++;
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
GlobalMatchingAnalyser::beginJob()
{
  edm::Service<TFileService> fs;
  TFileDirectory subDir = fs->mkdir( "matchAnalyzer" );

  //tg1 = new TGraphErrors();
  tg1 = subDir.make<TGraphErrors>(10); // TGraphErrors();
  tg1->SetName("tg1_name");
  tg1->SetTitle("tg1_title");
  tg1->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  tg1->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  glb_combined = subDir.make<TGraphErrors>(4);
  glb_combined->SetName("glb_combined");
  glb_combined->SetTitle("Global Combined Muon");
  glb_combined->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  glb_combined->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  glb_combined->SetLineColor(4);

  glb_inner = subDir.make<TGraphErrors>(4);
  glb_inner->SetName("glb_inner");
  glb_inner->SetTitle("Global Inner Muon");
  glb_inner->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  glb_inner->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  glb_inner->SetLineColor(3);

  glb_outer = subDir.make<TGraphErrors>(4);
  glb_outer->SetName("glb_outer");
  glb_outer->SetTitle("Global Outer Muon");
  glb_outer->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  glb_outer->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  glb_outer->SetLineColor(2);

  sta_muon = subDir.make<TGraphErrors>(4);
  sta_muon->SetName("sta_muon");
  sta_muon->SetTitle("Stand-alone Muon");
  sta_muon->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  sta_muon->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  sta_muon->SetLineColor(kRed);
  sta_muon->SetMarkerStyle(3);
  sta_muon->SetMarkerColor(2);

  tk_muon = subDir.make<TGraphErrors>(4);
  tk_muon->SetName("tk_muon");
  tk_muon->SetTitle("Tracker Muon");
  tk_muon->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  tk_muon->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  tk_muon->SetMarkerStyle(28);
  tk_muon->SetMarkerColor(3);

  region_fixed = subDir.make<TGraphErrors>(4);
  region_fixed->SetName("region_fixed");
  region_fixed->SetTitle("Fixed Region Size");
  region_fixed->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  region_fixed->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  region_dynamic = subDir.make<TGraphErrors>(4);
  region_dynamic->SetName("region_dynamic");
  region_dynamic->SetTitle("Dynamic Region Size");
  region_dynamic->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  region_dynamic->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  pos_tkAll = subDir.make<TGraphErrors>(10);
  pos_tkAll->SetName("pos_tkAll");
  pos_tkAll->SetTitle("All TK Tracks");
  pos_tkAll->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_tkAll->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  pos_tkAll->SetLineColor(kBlue);

  pos_tkCand = subDir.make<TGraphErrors>(4);
  pos_tkCand->SetName("pos_tkCand");
  pos_tkCand->SetTitle("TK Candidates");
  pos_tkCand->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_tkCand->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  pos_tkCandFixed = subDir.make<TGraphErrors>(4);
  pos_tkCandFixed->SetName("pos_tkCandFixed");
  pos_tkCandFixed->SetTitle("TK Candidates");
  pos_tkCandFixed->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_tkCandFixed->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  pos_selectedTkCand = subDir.make<TGraphErrors>(4);
  pos_selectedTkCand->SetName("pos_selectedTkCand");
  pos_selectedTkCand->SetTitle("Matched TK Candidates");
  pos_selectedTkCand->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_selectedTkCand->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  pos_selectedTkCand->SetMarkerStyle(24);

  pos_selectedTkCandFixed = subDir.make<TGraphErrors>(4);
  pos_selectedTkCandFixed->SetName("pos_selectedTkCandFixed");
  pos_selectedTkCandFixed->SetTitle("Matched TK Candidates");
  pos_selectedTkCandFixed->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_selectedTkCandFixed->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  pos_selectedTkCandFixed->SetMarkerStyle(24);
  
  surface_error1 = subDir.make<TGraphErrors>(4);
  surface_error1->SetName("surface_error1");
  surface_error1->SetTitle("Common Surface1");
  //surface_error1->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  //surface_error1->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  surface_error2 = subDir.make<TGraphErrors>(4);
  surface_error2->SetName("surface_error2");
  surface_error2->SetTitle("Common Surface2");
  //surface_error2->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  //surface_error2->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  
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


//define this as a plug-in
DEFINE_FWK_MODULE(GlobalMatchingAnalyser);
