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
// $Id: GlobalMatchingAnalyser.cc,v 1.5 2009/12/19 07:18:10 aeverett Exp $
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

  TFile *theFile; // self-explanatory
  TGraphErrors *tg1;
  TGraphErrors *glb_combined, *glb_inner, *glb_outer;
  TGraphErrors *sta_muon;
  TGraphErrors *tk_muon;

  TGraphErrors *region_dynamic, *region_fixed;

  TGraphErrors *pos_tkCand, *pos_tkCandFixed;
  TGraphErrors *pos_selectedTkCand, *pos_selectedTkCandFixed;

  TGraphErrors *surface_error1;
  TGraphErrors *surface_error2;

  std::string outputFileName;

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
  outputFileName = iConfig.getUntrackedParameter<std::string>("outputFileName");
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
  if(theFile) delete theFile;
  if(theRegionBuilder) delete theRegionBuilder;
  if(theService) delete theService;
  if(theTrackMatcher) delete theTrackMatcher;

  if(tg1) delete tg1;
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
      
      tkPreCandColl = chooseRegionalTrackerTracks(staTrack,tkTrackCands,iSta);
      tkPreCandCollFixed = chooseRegionalTrackerTracksFixed(staTrack,tkTrackCands,iSta);

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
      for(vector<TrackCand>::const_iterator iTk = tkPreCandCollFixed.begin();
	  iTk != tkPreCandCollFixed.end(); ++iTk) {
	//TrackRef tkRef(tkPreCandCollFixed,position);
	//position++;
	pos_tkCandFixed->SetPoint(iTkFixed,iTk->second->eta(),iTk->second->phi());
	pos_tkCandFixed->SetPointError(iTkFixed,iTk->second->etaError(),iTk->second->phiError());

	std::pair<TrajectoryStateOnSurface, TrajectoryStateOnSurface>
	  tsosPair = theTrackMatcher->convertToTSOSMuHit(staCand,*iTk);

	cout << "first " << tsosPair.first.isValid() << " second " << tsosPair.second.isValid() << endl;

	if(!tsosPair.first.isValid() || !tsosPair.second.isValid()) continue;

	// calculate matching variables
	double distance = theTrackMatcher->match_d(tsosPair.first,tsosPair.second);
	double chi2 = theTrackMatcher->match_Chi2(tsosPair.first,tsosPair.second);
	double loc_chi2 = theTrackMatcher->match_dist(tsosPair.first,tsosPair.second);
	double deltaR = theTrackMatcher->match_Rpos(tsosPair.first,tsosPair.second);
	cout << "eta " << iTk->second->eta() << " phi " << iTk->second->phi() << endl; 
	cout << "distance " << distance << endl;
	cout << "chi2 " << chi2 << endl;
	cout << "loc_chi2 " << loc_chi2 << endl;
	cout << "deltaR " << deltaR << endl;

	cout << "dR1 " << fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta()<1.5*0.2) << endl;
	cout << "dR2 " << (fabs(deltaPhi(tsosPair.second.globalPosition().phi(),tsosPair.first.globalPosition().phi())) < 0.2) << endl;

	cout << "dR1redo " << (fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta()) < 1.5 * 0.2) << " " << fabs(tsosPair.second.globalPosition().eta()-tsosPair.first.globalPosition().eta()) << endl;

	cout << "dR2redo " << fabs(deltaPhi(tsosPair.second.globalPosition().phi(),tsosPair.first.globalPosition().phi()) ) << endl;
	
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
	pos_selectedTkCandFixed->SetPoint(iSelTkFixed,iTk->second->eta(),iTk->second->phi());
	pos_selectedTkCandFixed->SetPointError(iSelTkFixed,iTk->second->etaError(),iTk->second->phiError());
	iSelTkFixed++;
      }
      
      iSta++;
    }

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
  theFile = new TFile(outputFileName.c_str(),"recreate");
  theFile->cd();
  tg1 = new TGraphErrors();
  tg1->SetName("tg1_name");
  tg1->SetTitle("tg1_title");
  tg1->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  tg1->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  glb_combined = new TGraphErrors();
  glb_combined->SetName("glb_combined");
  glb_combined->SetTitle("Global Combined Muon");
  glb_combined->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  glb_combined->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  glb_combined->SetLineColor(4);

  glb_inner = new TGraphErrors();
  glb_inner->SetName("glb_inner");
  glb_inner->SetTitle("Global Inner Muon");
  glb_inner->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  glb_inner->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  glb_inner->SetLineColor(3);

  glb_outer = new TGraphErrors();
  glb_outer->SetName("glb_outer");
  glb_outer->SetTitle("Global Outer Muon");
  glb_outer->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  glb_outer->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  glb_outer->SetLineColor(2);

  sta_muon = new TGraphErrors();
  sta_muon->SetName("sta_muon");
  sta_muon->SetTitle("Stand-alone Muon");
  sta_muon->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  sta_muon->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  sta_muon->SetLineColor(kRed);
  sta_muon->SetMarkerStyle(3);
  sta_muon->SetMarkerColor(2);

  tk_muon = new TGraphErrors();
  tk_muon->SetName("tk_muon");
  tk_muon->SetTitle("Tracker Muon");
  tk_muon->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  tk_muon->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  tk_muon->SetMarkerStyle(28);
  tk_muon->SetMarkerColor(3);

  region_fixed = new TGraphErrors();
  region_fixed->SetName("region_fixed");
  region_fixed->SetTitle("Fixed Region Size");
  region_fixed->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  region_fixed->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  region_dynamic = new TGraphErrors();
  region_dynamic->SetName("region_dynamic");
  region_dynamic->SetTitle("Dynamic Region Size");
  region_dynamic->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  region_dynamic->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  pos_tkCand = new TGraphErrors();
  pos_tkCand->SetName("pos_tkCand");
  pos_tkCand->SetTitle("TK Candidates");
  pos_tkCand->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_tkCand->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  pos_tkCandFixed = new TGraphErrors();
  pos_tkCandFixed->SetName("pos_tkCandFixed");
  pos_tkCandFixed->SetTitle("TK Candidates");
  pos_tkCandFixed->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_tkCandFixed->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  pos_selectedTkCand = new TGraphErrors();
  pos_selectedTkCand->SetName("pos_selectedTkCand");
  pos_selectedTkCand->SetTitle("Matched TK Candidates");
  pos_selectedTkCand->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_selectedTkCand->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  pos_selectedTkCand->SetMarkerStyle(24);

  pos_selectedTkCandFixed = new TGraphErrors();
  pos_selectedTkCandFixed->SetName("pos_selectedTkCandFixed");
  pos_selectedTkCandFixed->SetTitle("Matched TK Candidates");
  pos_selectedTkCandFixed->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_selectedTkCandFixed->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  pos_selectedTkCandFixed->SetMarkerStyle(24);
  
  surface_error1 = new TGraphErrors();
  surface_error1->SetName("surface_error1");
  surface_error1->SetTitle("Common Surface1");
  //surface_error1->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  //surface_error1->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  surface_error2 = new TGraphErrors();
  surface_error2->SetName("surface_error2");
  surface_error2->SetTitle("Common Surface2");
  //surface_error2->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  //surface_error2->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GlobalMatchingAnalyser::endJob() {
 theFile->cd();

 tg1->Write("",TObject::kOverwrite);

 glb_combined->Write("",TObject::kOverwrite);
 glb_inner->Write("",TObject::kOverwrite);
 glb_outer->Write("",TObject::kOverwrite);
 sta_muon->Write("",TObject::kOverwrite);
 tk_muon->Write("",TObject::kOverwrite);

 region_fixed->Write("",TObject::kOverwrite);
 region_dynamic->Write("",TObject::kOverwrite);

 pos_tkCand->Write("",TObject::kOverwrite);
 pos_tkCandFixed->Write("",TObject::kOverwrite);
 pos_selectedTkCand->Write("",TObject::kOverwrite);
 pos_selectedTkCandFixed->Write("",TObject::kOverwrite);
 
 surface_error1->Write("",TObject::kOverwrite);
 surface_error2->Write("",TObject::kOverwrite);
 
 theFile->Close();
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
