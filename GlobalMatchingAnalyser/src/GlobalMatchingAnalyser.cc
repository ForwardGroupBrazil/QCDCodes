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
// $Id: GlobalMatchingAnalyser.cc,v 1.2 2009/12/18 20:12:59 aeverett Exp $
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

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/GlobalTrackingTools/interface/MuonTrackingRegionBuilder.h"
#include "RecoTracker/TkTrackingRegions/interface/RectangularEtaPhiTrackingRegion.h"


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
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  RectangularEtaPhiTrackingRegion defineRegionOfInterest(const reco::TrackRef&) const;
  
  std::vector<reco::Track> chooseRegionalTrackerTracksFixed(const reco::TrackRef& staCand,
							    const edm::View<reco::Track>& tkTs,
							    int iSta);
  
  std::vector<reco::Track> chooseRegionalTrackerTracks(const reco::TrackRef& staCand,
						       const edm::View<reco::Track>& tkTs,
						       int iSta);
  
  // ----------member data ---------------------------
  MuonTrackingRegionBuilder* theRegionBuilder;

  MuonServiceProxy* theService;

  TFile *theFile; // self-explanatory
  TGraphErrors *tg1;
  TGraphErrors *glb_combined, *glb_inner, *glb_outer;
  TGraphErrors *sta_muon;
  TGraphErrors *tk_muon;

  TGraphErrors *region_dynamic, *region_fixed;

  TGraphErrors *tkCand, *tkCandFixed;

  std::string outputFileName;

  edm::InputTag theTrackLabel;
  edm::InputTag theMuonLabel;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//
GlobalMatchingAnalyser::GlobalMatchingAnalyser(const edm::ParameterSet& iConfig) : theRegionBuilder(0), theService(0)

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


}


GlobalMatchingAnalyser::~GlobalMatchingAnalyser()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  if(theFile) delete theFile;
  if(theRegionBuilder) delete theRegionBuilder;
  if(theService) delete theService;

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

  tg1->SetPoint(0,1,1);
  tg1->SetPoint(1,2,1);
  tg1->SetPointError(0,1,.5);
  tg1->SetPointError(1,0.25,1.);

  int nMu = muonColl.size();

  int iMu = 0;
  int iSta = 0;
  int iTkFixed = 0;
  int iTkDynamic = 0;
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

    vector<Track> tkCandColl;
    vector<Track> tkCandCollFixed;

    if(iMuon->isStandAloneMuon()) {
      if(staTrack.isAvailable()) {
	sta_muon->SetPoint(iMu,staTrack->eta(),staTrack->phi());
	sta_muon->SetPointError(iMu,staTrack->etaError(),staTrack->phiError());
      }
      
      tkCandColl = chooseRegionalTrackerTracks(staTrack,trkColl,iSta);
      tkCandCollFixed = chooseRegionalTrackerTracksFixed(staTrack,trkColl,iSta);

      for(vector<Track>::const_iterator iTk = tkCandColl.begin();
	  iTk != tkCandColl.end(); ++iTk) {
	tkCand->SetPoint(iTkDynamic,iTk->eta(),iTk->phi());
	tkCand->SetPointError(iTkDynamic,iTk->etaError(),iTk->phiError());
	iTkDynamic++;
      }

      for(vector<Track>::const_iterator iTk = tkCandCollFixed.begin();
	  iTk != tkCandCollFixed.end(); ++iTk) {
	tkCandFixed->SetPoint(iTkFixed,iTk->eta(),iTk->phi());
	tkCandFixed->SetPointError(iTkFixed,iTk->etaError(),iTk->phiError());
	iTkFixed++;
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

  glb_inner = new TGraphErrors();
  glb_inner->SetName("glb_inner");
  glb_inner->SetTitle("Global Inner Muon");
  glb_inner->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  glb_inner->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  glb_outer = new TGraphErrors();
  glb_outer->SetName("glb_outer");
  glb_outer->SetTitle("Global Outer Muon");
  glb_outer->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  glb_outer->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  sta_muon = new TGraphErrors();
  sta_muon->SetName("sta_muon");
  sta_muon->SetTitle("Stand-alone Muon");
  sta_muon->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  sta_muon->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  tk_muon = new TGraphErrors();
  tk_muon->SetName("tk_muon");
  tk_muon->SetTitle("Tracker Muon");
  tk_muon->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  tk_muon->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

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

  tkCand = new TGraphErrors();
  tkCand->SetName("tkCand");
  tkCand->SetTitle("TK Candidates");
  tkCand->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  tkCand->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

  tkCandFixed = new TGraphErrors();
  tkCandFixed->SetName("tkCandFixed");
  tkCandFixed->SetTitle("TK Candidates");
  tkCandFixed->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  tkCandFixed->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);

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

 tkCand->Write("",TObject::kOverwrite);
 tkCandFixed->Write("",TObject::kOverwrite);

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
vector<Track> 
GlobalMatchingAnalyser::chooseRegionalTrackerTracksFixed(const reco::TrackRef& staCand, 
                                                         const View<Track>& tkTs,
							 int iSta) {
  
  // define eta-phi region
  RectangularEtaPhiTrackingRegion regionOfInterest = defineRegionOfInterest(staCand);
  
  // get region's etaRange and phiMargin
  PixelRecoRange<float> etaRange = regionOfInterest.etaRange();
  TkTrackingRegionsMargin<float> phiMargin = regionOfInterest.phiMargin();

  region_fixed->SetPoint(iSta,regionOfInterest.direction().eta(),regionOfInterest.direction().phi());
  region_fixed->SetPointError(iSta,1.,1.);

  vector<Track> result;

  double deltaR_max = 1.0;

  for ( View<Track>::const_iterator is = tkTs.begin(); 
	is != tkTs.end(); ++is ) {

    double deltaR_tmp = deltaR(static_cast<double>(regionOfInterest.direction().eta()),
			       static_cast<double>(regionOfInterest.direction().phi()),
			       is->eta(), is->phi());

    // for each trackCand in region, add trajectory and add to result
    if (deltaR_tmp < deltaR_max) {
      //Track tmpCand = TrackCand(*is);
      //result.push_back(tmpCand);
      result.push_back(*is);
    }
  }

  return result; 

}

//
// select tracker tracks within a region of interest
//
vector<Track> 
GlobalMatchingAnalyser::chooseRegionalTrackerTracks( const TrackRef& staCand, 
						     const View<Track>& tkTs,
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

  vector<Track> result;

  for ( View<Track>::const_iterator is = tkTs.begin(); 
	is != tkTs.end(); ++is ) {
    
    double trackEta = is->eta();
    double trackPhi = is->phi();
    
    // Clean  
    bool inEtaRange = trackEta >= (etaRange.mean() - etaLimit) && trackEta <= (etaRange.mean() + etaLimit) ;
    bool inPhiRange = (fabs(deltaPhi(trackPhi,double(region.direction().phi()))) < phiLimit );
    
    if(inEtaRange && inPhiRange) result.push_back(*is);
    
  }

  return result;

}


//define this as a plug-in
DEFINE_FWK_MODULE(GlobalMatchingAnalyser);
