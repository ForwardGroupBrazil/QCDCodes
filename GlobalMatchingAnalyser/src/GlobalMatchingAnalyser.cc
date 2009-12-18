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
// $Id$
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

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include <DataFormats/TrackReco/interface/Track.h>

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
  
  // ----------member data ---------------------------
  TFile *theFile; // self-explanatory
  TGraphErrors *tg1;
  TGraphErrors *glb_combined, *glb_inner, *glb_outer;
  TGraphErrors *sta_muon;
  TGraphErrors *tk_muon;

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

//
// constructors and destructor
//
GlobalMatchingAnalyser::GlobalMatchingAnalyser(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  outputFileName = iConfig.getUntrackedParameter<std::string>("outputFileName");
  theTrackLabel = iConfig.getParameter<edm::InputTag>("trackLabel");
  theMuonLabel = iConfig.getParameter<edm::InputTag>("muonLabel");

}


GlobalMatchingAnalyser::~GlobalMatchingAnalyser()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  if(theFile) delete theFile;
  if(tg1) delete tg1;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GlobalMatchingAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;



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
  for(View<Muon>::const_iterator iMuon = muonColl.begin();
      iMuon != muonColl.end(); ++iMuon) {

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

    if(iMuon->isStandAloneMuon()) {
      if(staTrack.isAvailable()) {
	sta_muon->SetPoint(iMu,staTrack->eta(),staTrack->phi());
	sta_muon->SetPointError(iMu,staTrack->etaError(),staTrack->phiError());
      }
    }

    if(iMuon->isTrackerMuon()) {
      if(tkTrack.isAvailable()) {
	tk_muon->SetPoint(iMu,tkTrack->eta(),tkTrack->phi());
	tk_muon->SetPointError(iMu,tkTrack->etaError(),tkTrack->phiError());
      }
    }

    iMu++;
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
 theFile->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(GlobalMatchingAnalyser);
