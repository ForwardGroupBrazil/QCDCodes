// -*- C++ -*-
//
// Package:    RecoMuonEffAnlyzer
// Class:      RecoMuonEffAnlyzer
// 
/**\class RecoMuonEffAnlyzer RecoMuonEffAnlyzer.cc AEverett/RecoMuonEffAnlyzer/src/RecoMuonEffAnlyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam A Everett
//         Created:  Wed Sep 27 14:54:28 EDT 2006
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
//#include <TH2.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFrame.h>
//#include <TMath.h>
//#include <TF1.h>

//
// class decleration
//

using namespace std;
using namespace edm;

class RecoMuonEffAnlyzer : public edm::EDAnalyzer {
public:
  explicit RecoMuonEffAnlyzer(const edm::ParameterSet&);
  ~RecoMuonEffAnlyzer();
  
  typedef std::pair<reco::Track, reco::Track> TrackCand;
  typedef std::pair<reco::Track, SimTrack> TrackCandSim;
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual float calculateDistance(const math::XYZVector&, const math::XYZVector&);
  virtual TH1F* divideErr(TH1F*, TH1F*, TH1F*);

  // ----------member data ---------------------------
  string out, open;
  double  theMinEta, theMaxEta, theMinPt;
  int theNBins, thePartID;

  edm::InputTag GLBtrackTags_; 
  edm::InputTag STAtrackTags_; 
  edm::InputTag MuonTags_; 
  edm::InputTag SIMtrackTags_; 

  Handle<reco::MuonCollection> MuCollection;
  Handle<reco::TrackCollection> GLBTrackCollection;
  Handle<reco::TrackCollection> STATrackCollection;
  Handle<edm::SimTrackContainer> SIMTrackCollection;
  
  MuonServiceProxy* theService;

  //ROOT Pointers
  TFile* hFile;
  TStyle* effStyle;

  TH1F* hi_dist_glb_sta;
  TH1F* hi_dist_glb_sim;
  TH1F* hi_dist_sta_sim;

  TH1F* hi_sim_pt  ;
  TH1F* hi_sta_pt  ;
  TH1F* hi_glb_pt  ;

  TH1F* hi_sim_eta  ;
  TH1F* hi_sta_eta  ;
  TH1F* hi_glb_eta  ;

  TH1F* hi_glbsim_pt  ;
  TH1F* hi_stasim_pt  ;
  TH1F* hi_glbsta_pt  ;

  TH1F* hi_glbsim_eta  ;
  TH1F* hi_stasim_eta  ;
  TH1F* hi_glbsta_eta  ;

  TH1F* hi_glbsim2_pt  ;
  TH1F* hi_stasim2_pt  ;
  TH1F* hi_glbsta2_pt  ;

  TH1F* hi_glbsim2_eta  ;
  TH1F* hi_stasim2_eta  ;
  TH1F* hi_glbsta2_eta  ;

  TH1F* hi_glbsim_eff_pt  ;
  TH1F* hi_stasim_eff_pt  ;
  TH1F* hi_glbsta_eff_pt  ;

  TH1F* hi_glbsim_eff_eta  ;
  TH1F* hi_stasim_eff_eta  ;
  TH1F* hi_glbsta_eff_eta  ;

  TH1F* hi_glbsim_pur_pt  ;
  TH1F* hi_stasim_pur_pt  ;
  TH1F* hi_glbsta_pur_pt  ;

  TH1F* hi_glbsim_pur_eta  ;
  TH1F* hi_stasim_pur_eta  ;
  TH1F* hi_glbsta_pur_eta  ;

  TH1F* hi_glbsim_ptres;
  TH1F* hi_stasim_ptres;
  TH1F* hi_glbsta_ptres;

  TH1F* hi_glbsim_etares;
  TH1F* hi_stasim_etares;
  TH1F* hi_glbsta_etares;
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
RecoMuonEffAnlyzer::RecoMuonEffAnlyzer(const edm::ParameterSet& iConfig) 
  :
  GLBtrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("GLBtracks")),
  STAtrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("STAtracks")),
  MuonTags_(iConfig.getUntrackedParameter<edm::InputTag>("Muons")),
  SIMtrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("SIMtracks")),
  out(iConfig.getParameter<string>("out")),
  open(iConfig.getParameter<string>("open")),
  theMinEta(iConfig.getParameter<double>("etamin")),
  theMaxEta(iConfig.getParameter<double>("etamax")),
  theMinPt(iConfig.getParameter<double>("ptmin")),
  theNBins(iConfig.getParameter<int>("nbins")),
  thePartID(iConfig.getParameter<int>("partId"))
{
  //now do what ever initialization is needed

  // service parameters
  ParameterSet serviceParameters = iConfig.getParameter<ParameterSet>("ServiceParameters");
  // the services
  theService = new MuonServiceProxy(serviceParameters);
  
}


RecoMuonEffAnlyzer::~RecoMuonEffAnlyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  if (hFile!=0) {
    hFile->Close();
    delete hFile;
  }
  if (theService) delete theService;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
RecoMuonEffAnlyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  cout << "DEBUG  210" << endl;
  //using namespace edm;
  using reco::TrackCollection;
  using reco::MuonCollection;
  cout << "DEBUG  214" << endl;
  // Update the services
  //theService->update(iSetup);
  
  //   Handle<TrackCollection> GLBTrackCollection;
  iEvent.getByLabel( GLBtrackTags_, GLBTrackCollection);
  const reco::TrackCollection glbTC = *(GLBTrackCollection.product());
  cout << "GLB size: " << glbTC.size() << endl;
  //   Handle<TrackCollection> STATrackCollection;
  iEvent.getByLabel( STAtrackTags_, STATrackCollection);
  const reco::TrackCollection staTC = *(STATrackCollection.product());
  cout << "STA size: " << staTC.size() << endl;cout << "DEBUG  219" << endl;
  //Handle<reco::MuonCollection> MuCollection;cout << "DEBUG  220" << endl;
  iEvent.getByLabel(MuonTags_,MuCollection);cout << "DEBUG  221" << endl;
  const reco::MuonCollection muonC = *(MuCollection.product());cout << "DEBUG  222" << endl;
  cout << "Mu size: " << muonC.size() << endl;
  //   Handle<SimTrackContainer> SIMTrackCollection;
  iEvent.getByLabel(SIMtrackTags_,SIMTrackCollection);cout << "DEBUG  231" << endl;
  const SimTrackContainer simTC = *(SIMTrackCollection.product());
  cout << "DEBUG  227" << endl;
  vector<SimTrack> simMuons1, simMuons2;
  //int position = 0;
  SimTrackContainer::const_iterator simTrack;
  for (simTrack = simTC.begin(); simTrack != simTC.end(); ++simTrack){
    if (abs((*simTrack).type()) == 13
	&& (*simTrack).momentum().perp() >= theMinPt 
	&& (*simTrack).momentum().eta() <= theMaxEta
	&& (*simTrack).momentum().eta() >= theMinEta
	) {
      //position++;
      cout<<"Sim pT: "<<(*simTrack).momentum().perp()<<endl;
      //simPt=(*simTrack).momentum().perp();
      cout<<"Sim Eta: "<<(*simTrack).momentum().eta()<<endl;
      //numberOfSimTracks++;
      //SimTrackRef simTrackRef(simTC,position-1);
      simMuons1.push_back(*simTrack);
      simMuons2.push_back(*simTrack);
      hi_sim_pt->Fill((*simTrack).momentum().perp());
      hi_sim_eta->Fill(((*simTrack).momentum().eta()));
    }    
  }
  cout << "DEBUG  249" << endl;
  vector<reco::Track> staMuons1, staMuons2;
  TrackCollection::const_iterator staTrack;
  for (staTrack = staTC.begin(); staTrack != staTC.end(); ++staTrack){
    cout<<"STA pT: "<<(*staTrack).pt()<<endl;
    //simPt=(*simTrack).momentum().perp();
    cout<<"STA Eta: "<<(*staTrack).eta()<<endl;
    staMuons1.push_back(*staTrack);
    staMuons2.push_back(*staTrack);
    hi_sta_pt->Fill((*staTrack).pt());
    hi_sta_eta->Fill(((*staTrack).eta()));         
  }
  cout << "DEBUG  261" << endl;
  vector<reco::Track> glbMuons1, glbMuons2;
  TrackCollection::const_iterator glbTrack;
  for (glbTrack = glbTC.begin(); glbTrack != glbTC.end(); ++glbTrack){
    cout<<"GLB pT: "<<(*glbTrack).pt()<<endl;
    //simPt=(*simTrack).momentum().perp();
    cout<<"GLB Eta: "<<(*glbTrack).eta()<<endl;
    glbMuons1.push_back(*glbTrack);
    glbMuons2.push_back(*glbTrack);
    hi_glb_pt->Fill((*glbTrack).pt());
    hi_glb_eta->Fill(((*glbTrack).eta()));         
  }
  cout << "DEBUG  273" << endl;
  for (glbTrack = glbTC.begin(); glbTrack != glbTC.end(); ++glbTrack){
    for (staTrack = staTC.begin(); staTrack != staTC.end(); ++staTrack){
      float D = calculateDistance(glbTrack->momentum(),staTrack->momentum());
      hi_dist_glb_sta->Fill(D);
    }
  }
  
  for (glbTrack = glbTC.begin(); glbTrack != glbTC.end(); ++glbTrack){
    for (vector<SimTrack>::const_iterator simTrack = simMuons1.begin(); simTrack != simMuons1.end(); ++simTrack){
      const math::XYZVector simVect(simTrack->momentum().x(),simTrack->momentum().y(),simTrack->momentum().z());
      float D = calculateDistance(glbTrack->momentum(),simVect);
      hi_dist_glb_sim->Fill(D);
    }
  }
  
  for (staTrack = staTC.begin(); staTrack != staTC.end(); ++staTrack){
    for (vector<SimTrack>::const_iterator simTrack = simMuons1.begin(); simTrack != simMuons1.end(); ++simTrack){
      const math::XYZVector simVect(simTrack->momentum().x(),simTrack->momentum().y(),simTrack->momentum().z());
      float D = calculateDistance(staTrack->momentum(),simVect);
      hi_dist_sta_sim->Fill(D);
    }
  }
  
  //Now make pairs
  //FIXME: add a Distance cut
  
  //Global-Simulated pairs
  vector<TrackCandSim> pairGlbSim;
  for (vector<SimTrack>::const_iterator simTrack = simMuons1.begin(); simTrack != simMuons1.end(); ++simTrack){
    float minDist = 1000.;
    int index = 0;
    int keep = -1;
    cout << "glbMuons1.size = " << glbMuons1.size() << endl;
    for (glbTrack = glbMuons1.begin(); glbTrack != glbMuons1.end(); ++glbTrack){
      const math::XYZVector simVect(simTrack->momentum().x(),simTrack->momentum().y(),simTrack->momentum().z());
      float D = calculateDistance(glbTrack->momentum(),simVect); 
      if(D <= minDist) {
	minDist = D;
	keep = index;
      }
      index++;
    }
    if(keep > -1 ) {
      reco::Track tmp(glbMuons1.at(keep));
      TrackCandSim tkCandSim = TrackCandSim(tmp,*simTrack);
      pairGlbSim.push_back(tkCandSim);
      cout << "pairGlbSim.size = " << pairGlbSim.size() << " index " << keep << endl;
      glbMuons1.erase(glbMuons1.begin()+keep,glbMuons1.begin()+keep+1);
    }
  }
  
  for (vector<TrackCandSim>::const_iterator pair = pairGlbSim.begin(); pair != pairGlbSim.end(); pair++){
    hi_glbsim_pt->Fill((*pair).second.momentum().perp());
    hi_glbsim_eta->Fill(((*pair).second.momentum().eta()));       
    hi_glbsim2_pt->Fill((*pair).first.pt());
    hi_glbsim2_eta->Fill(((*pair).first.eta()));       
    
    float ptSim = (*pair).second.momentum().perp();
    float ptGLB = (*pair).first.pt();
    hi_glbsim_ptres->Fill((1/ptGLB - 1/ptSim)/(1/ptSim));
    
    float etaSim = (*pair).second.momentum().eta();
    float etaGLB = (*pair).first.eta();
    hi_glbsim_etares->Fill((1/etaGLB - 1/etaSim)/(1/etaSim));
  }
  
  //StandAlone-Simulated pair
  vector<TrackCandSim> pairStaSim;
  for (vector<SimTrack>::const_iterator simTrack = simMuons1.begin(); simTrack != simMuons1.end(); ++simTrack){
    float minDist = 1000.;
    int index = 0;
    int keep = -1;
    cout << "staMuons1.size = " << staMuons1.size() << endl;
    for (staTrack = staMuons1.begin(); staTrack != staMuons1.end(); ++staTrack){
      const math::XYZVector simVect(simTrack->momentum().x(),simTrack->momentum().y(),simTrack->momentum().z());
      float D = calculateDistance(staTrack->momentum(),simVect); 
      if(D <= minDist) {
	minDist = D;
	keep = index;
      }
      index++;
    }
    if(keep > -1 ) {
      reco::Track tmp(staMuons1.at(keep));
      TrackCandSim tkCandSim = TrackCandSim(tmp,*simTrack);
      pairStaSim.push_back(tkCandSim);
      cout << "pairStaSim.size = " << pairStaSim.size() << " index " << keep << endl;
      staMuons1.erase(staMuons1.begin()+keep,staMuons1.begin()+keep+1);
    }cout << "Debug 362" << endl;
  }
  cout << "Debug 364" << endl;
  for (vector<TrackCandSim>::const_iterator pair = pairStaSim.begin(); pair != pairStaSim.end(); pair++){cout << "Debug 365" << endl;
  hi_stasim_pt->Fill((*pair).second.momentum().perp());cout << "Debug 366" << endl;
  hi_stasim_eta->Fill(((*pair).second.momentum().eta()));cout << "Debug 367" << endl;       
  hi_stasim2_pt->Fill((*pair).first.pt());
  hi_stasim2_eta->Fill(((*pair).first.eta()));
  
  float ptSim = (*pair).second.momentum().perp();cout << "Debug 371" << endl;
  float ptSTA = (*pair).first.pt();cout << "Debug 372" << endl;
  hi_stasim_ptres->Fill((1/ptSTA-1/ptSim)/(1/ptSim));cout << "Debug 372" << endl;
  
  float etaSim = (*pair).second.momentum().eta();cout << "Debug 375" << endl;
  float etaSTA = (*pair).first.eta();cout << "Debug 364" << endl;
  hi_stasim_etares->Fill((1/etaSTA-1/etaSim)/(1/etaSim));cout << "Debug 377" << endl;
  }
  cout << "Debug 386" << endl;
  //Global-StandAlone pairs
  
  for (reco::MuonCollection::const_iterator pair = muonC.begin(); pair != muonC.end(); pair++){ 
    reco::TrackRef global = pair->combinedMuon();
    reco::TrackRef standAlone = pair->standAloneMuon();
    standAlone->pt();
    //pair->combinedMuon()->momentum().perp2();
    hi_glbsta_pt->Fill( sqrt(global->pt()) );
    hi_glbsta_eta->Fill( (global->eta() ) );
    //hi_glbsta_pt->Fill( sqrt(standAlone->pt()) );
    //hi_glbsta_eta->Fill( (standAlone->eta() ) );
  }
  
  
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
RecoMuonEffAnlyzer::beginJob(const edm::EventSetup&)
{
   hFile = new TFile( out.c_str(), open.c_str() );

   effStyle = new TStyle("effStyle","Efficiency Study Style");
   
   //gROOT->Reset();cout << "DEBUG  404" << endl;
   effStyle->SetCanvasBorderMode(0);cout << "DEBUG  405" << endl;
   //   effStyle->SetCanvasBorderMode(0);cout << "DEBUG  406" << endl;
   effStyle->SetOptTitle(0);cout << "DEBUG  407" << endl;
   effStyle->SetOptStat(0);cout << "DEBUG  408" << endl;
   //...................  
   effStyle->SetCanvasBorderMode(0);cout << "Debug 410" << endl;
   effStyle->SetPadBorderMode(1);cout << "Debug 411" << endl;
   effStyle->SetOptTitle(0);cout << "Debug 412" << endl;
   effStyle->SetStatFont(42);cout << "Debug 413" << endl;
   effStyle->SetTitleFont(22);cout << "Debug 414" << endl;
   effStyle->SetCanvasColor(10);cout << "Debug 415" << endl;
   effStyle->SetPadColor(0);cout << "Debug 416" << endl;
   effStyle->SetLabelFont(42,"x");cout << "Debug 417" << endl;
   effStyle->SetLabelFont(42,"y");cout << "Debug 418" << endl;
   effStyle->SetHistFillStyle(1001);cout << "Debug 419" << endl;
   effStyle->SetHistFillColor(0);cout << "Debug 420" << endl;
   effStyle->SetOptStat(0);cout << "Debug 421" << endl;
   //effStyle->SetOptStat(00011110);cout << "Debug 422" << endl;
   effStyle->SetOptFit(0111);cout << "Debug 423" << endl;
   effStyle->SetStatH(0.05); cout << "Debug 424" << endl;
   //.................... 
   

   hi_dist_glb_sta = new TH1F("hi_dist_glb_sta","Distance GLB-STA",100,0.,100.);
   hi_dist_glb_sim = new TH1F("hi_dist_glb_sim","Distance GLB-SIM",100,0.,100.);
   hi_dist_sta_sim = new TH1F("hi_dist_sta_sim","Distance STA-SIM",100,0.,100.);

   hi_sim_pt  = new TH1F("hi_sim_pt","P_{T}^{sim}",theNBins,0.0,100.);
   hi_sim_pt->Sumw2();
   hi_sta_pt   = new TH1F("hi_sta_pt","P_{T}^{STA}",theNBins,0.0,100.);
   hi_sta_pt->Sumw2();
   hi_glb_pt   = new TH1F("hi_glb_pt","P_{T}^{GLB}",theNBins,0.0,100.);
   hi_glb_pt->Sumw2();

   hi_sim_eta  = new TH1F("hi_sim_eta","#eta^{sim}",theNBins,theMinEta,theMaxEta);
   hi_sim_eta->Sumw2();
   hi_sta_eta   = new TH1F("hi_sta_eta","#eta^{STA}",theNBins,theMinEta,theMaxEta);
   hi_sta_eta->Sumw2();
   hi_glb_eta   = new TH1F("hi_glb_eta","#eta^{GLB}",theNBins,theMinEta,theMaxEta);
   hi_glb_eta->Sumw2();

   hi_glbsim_pt  = new TH1F("hi_glbsim_pt","P_{T}^{GLB,sim}",theNBins,0.0,100.);
   hi_glbsim_pt->Sumw2();
   hi_stasim_pt   = new TH1F("hi_stasim_pt","P_{T}^{STA,sim}",theNBins,0.0,100.);
   hi_stasim_pt->Sumw2();
   hi_glbsta_pt   = new TH1F("hi_glbsta_pt","P_{T}^{GLB,STA}",theNBins,0.0,100.);
   hi_glbsta_pt->Sumw2();

   hi_glbsim_eta  = new TH1F("hi_glbsim_eta","#eta^{GLB,sim}",theNBins,theMinEta,theMaxEta);
   hi_glbsim_eta->Sumw2();
   hi_stasim_eta   = new TH1F("hi_stasim_eta","#eta^{STA,sim}",theNBins,theMinEta,theMaxEta);
   hi_stasim_eta->Sumw2();
   hi_glbsta_eta   = new TH1F("hi_glbsta_eta","#eta^{GLB,STA}",theNBins,theMinEta,theMaxEta);
   hi_glbsta_eta->Sumw2();

   hi_glbsim2_pt  = new TH1F("hi_glbsim2_pt","P_{T}^{GLB,sim}",theNBins,0.0,100.);
   hi_glbsim2_pt->Sumw2();
   hi_stasim2_pt   = new TH1F("hi_stasim2_pt","P_{T}^{STA,sim}",theNBins,0.0,100.);
   hi_stasim2_pt->Sumw2();
   hi_glbsta2_pt   = new TH1F("hi_glbsta2_pt","P_{T}^{GLB,STA}",theNBins,0.0,100.);
   hi_glbsta2_pt->Sumw2();

   hi_glbsim2_eta  = new TH1F("hi_glbsim2_eta","#eta^{GLB,sim}",theNBins,theMinEta,theMaxEta);
   hi_glbsim2_eta->Sumw2();
   hi_stasim2_eta   = new TH1F("hi_stasim2_eta","#eta^{STA,sim}",theNBins,theMinEta,theMaxEta);
   hi_stasim2_eta->Sumw2();
   hi_glbsta2_eta   = new TH1F("hi_glbsta2_eta","#eta^{GLB,STA}",theNBins,theMinEta,theMaxEta);
   hi_glbsta2_eta->Sumw2();

   hi_glbsim_eff_pt  = new TH1F("hi_glbsim_eff_pt","Efficiency GLB,sim",theNBins,0.0,100.);
   hi_stasim_eff_pt   = new TH1F("hi_stasim_eff_pt","Efficiency STA,sim",theNBins,0.0,100.);
   hi_glbsta_eff_pt   = new TH1F("hi_glbsta_eff_pt","Efficiency GLB,STA",theNBins,0.0,100.);

   hi_glbsim_eff_eta  = new TH1F("hi_glbsim_eff_eta","Efficiency GLB,sim",theNBins,theMinEta,theMaxEta);
   hi_stasim_eff_eta   = new TH1F("hi_stasim_eff_eta","Efficiency STA,sim",theNBins,theMinEta,theMaxEta);
   hi_glbsta_eff_eta   = new TH1F("hi_glbsta_eff_eta","Efficiency GLB,STA",theNBins,theMinEta,theMaxEta);

   hi_glbsim_pur_pt  = new TH1F("hi_glbsim_pur_pt","Purity GLB,sim",theNBins,0.0,100.);
   hi_stasim_pur_pt   = new TH1F("hi_stasim_pur_pt","Purity STA,sim",theNBins,0.0,100.);
   hi_glbsta_pur_pt   = new TH1F("hi_glbsta_pur_pt","Purity GLB,STA",theNBins,0.0,100.);

   hi_glbsim_pur_eta  = new TH1F("hi_glbsim_pur_eta","Purity GLB,sim",theNBins,theMinEta,theMaxEta);
   hi_stasim_pur_eta   = new TH1F("hi_stasim_pur_eta","Purity STA,sim",theNBins,theMinEta,theMaxEta);
   hi_glbsta_pur_eta   = new TH1F("hi_glbsta_pur_eta","Purity GLB,STA",theNBins,theMinEta,theMaxEta);


  hi_glbsim_ptres  = new TH1F("hi_glbsim_ptres","P_{T} resolution GLB,sim",theNBins,-0.2,0.2);
  hi_stasim_ptres  = new TH1F("hi_stasim_ptres","P_{T} resolution STA.sim",theNBins,-0.2,0.2);
  hi_glbsta_ptres  = new TH1F("hi_glbsta_ptres","P_{T} resolution GLB,STA",theNBins,-0.2,0.2);

  hi_glbsim_etares  = new TH1F("hi_glbsim_etares","#eta resolution GLB,sim",theNBins,-1.,1.);
  hi_stasim_etares  = new TH1F("hi_stasim_etares","#eta resolution STA.sim",theNBins,-1.,1.);
  hi_glbsta_etares  = new TH1F("hi_glbsta_etares","#eta resolution GLB,STA",theNBins,-0.1,0.1);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RecoMuonEffAnlyzer::endJob() {
  cout << "Debug 503" << endl;
  hi_glbsim_eff_pt = divideErr(hi_glbsim_pt,hi_sim_pt,hi_glbsim_eff_pt);
  hi_stasim_eff_pt = divideErr(hi_stasim_pt,hi_sim_pt,hi_stasim_eff_pt);
  hi_glbsta_eff_pt = divideErr(hi_glbsta_pt,hi_sta_pt,hi_glbsta_eff_pt);

  hi_glbsim_eff_eta = divideErr(hi_glbsim_eta,hi_sim_eta,hi_glbsim_eff_eta);
  hi_stasim_eff_eta = divideErr(hi_stasim_eta,hi_sim_eta,hi_stasim_eff_eta);
  hi_glbsta_eff_eta = divideErr(hi_glbsta_eta,hi_sta_eta,hi_glbsta_eff_eta);

  //FIXME: there is a problem with migration
  cout << "Debug 513" << endl;
  hi_glbsim_pur_pt = divideErr(hi_glbsim_pt,hi_glb_pt,hi_glbsim_pur_pt);
  hi_stasim_pur_pt = divideErr(hi_stasim_pt,hi_sta_pt,hi_stasim_pur_pt);
  hi_glbsta_pur_pt = divideErr(hi_glbsta_pt,hi_glb_pt,hi_glbsta_pur_pt);
  cout << "Debug 517" << endl;
  hi_glbsim_pur_eta = divideErr(hi_glbsim_eta,hi_glb_eta,hi_glbsim_pur_eta);
  hi_stasim_pur_eta = divideErr(hi_stasim_eta,hi_sta_eta,hi_stasim_pur_eta);
  hi_glbsta_pur_eta = divideErr(hi_glbsta_eta,hi_glb_eta,hi_glbsta_pur_eta);
  cout << "Debug 521" << endl;
  //...................  
  effStyle->SetCanvasBorderMode(0);cout << "Debug 523" << endl;
  effStyle->SetPadBorderMode(1);cout << "Debug 524" << endl;
  effStyle->SetOptTitle(0);cout << "Debug 525" << endl;
  effStyle->SetStatFont(42);cout << "Debug 526" << endl;
  effStyle->SetTitleFont(22);cout << "Debug 527" << endl;
  effStyle->SetCanvasColor(10);cout << "Debug 528" << endl;
  effStyle->SetPadColor(0);cout << "Debug 529" << endl;
  effStyle->SetLabelFont(42,"x");cout << "Debug 530" << endl;
  effStyle->SetLabelFont(42,"y");cout << "Debug 531" << endl;
  effStyle->SetHistFillStyle(1001);cout << "Debug 532" << endl;
  effStyle->SetHistFillColor(0);cout << "Debug 533" << endl;
  effStyle->SetOptStat(0);cout << "Debug 534" << endl;
  //effStyle->SetOptStat(00011110);cout << "Debug 535" << endl;
  effStyle->SetOptFit(0111);cout << "Debug 536" << endl;
  effStyle->SetStatH(0.05); cout << "Debug 537" << endl;
  //.................... 
  cout << "Debug 539" << endl;
  gROOT->SetStyle("effStyle");
  cout << "Debug 541" << endl;
  char* l1string = "Simulated";
  char* l2string = "StandAlone";
  char* l3string = "Global";
  cout << "Debug 545" << endl;

  TCanvas* c1 = new TCanvas("eff","Efficiency pt",10,10,700,500);
  c1->SetFillColor(0);
  c1->SetGrid(1);
  c1->SetTicky();
  c1->SetRightMargin(0.03);
  c1->SetTopMargin(0.02);
  c1->cd(); 
  cout << "Debug 554" << endl;  
  hi_stasim_eff_pt->SetXTitle("p_{T}^{#mu}");
  hi_stasim_eff_pt->SetYTitle("Efficiency");
  hi_stasim_eff_pt->SetTitleOffset(1.1,"x");
  hi_stasim_eff_pt->SetTitleOffset(1.15,"y");
  //hi_stasim_eff_pt->SetMaximum(1.02);
  //hi_stasim_eff_pt->SetMinimum(0.5);
  hi_stasim_eff_pt->SetStats(false);
  cout << "Debug 562" << endl;
  hi_stasim_eff_pt->SetLineWidth(1);
  hi_glbsim_eff_pt->SetLineWidth(1);
  hi_stasim_eff_pt->SetLineColor(2);
  hi_glbsim_eff_pt->SetLineColor(4);
  hi_stasim_eff_pt->SetLineStyle(1);
  hi_glbsim_eff_pt->SetLineStyle(1);

  //hi_stasim_eff_pt->SetFillColor(5);
  //hi_glbsim_eff_pt->SetFillColor(7);

  hi_stasim_eff_pt->SetMarkerStyle(21);
  hi_glbsim_eff_pt->SetMarkerStyle(26);
  hi_stasim_eff_pt->SetMarkerColor(2);
  hi_glbsim_eff_pt->SetMarkerColor(4);
  cout << "Debug 577" << endl;
  hi_stasim_eff_pt ->DrawCopy("PE");
  hi_glbsim_eff_pt ->DrawCopy("PEsame");
  //hi_stasim_eff_pt ->DrawCopy("AxisSame");
  cout << "Debug 581" << endl;
  TLegend* legend1 = new TLegend(0.6,0.2,0.8,0.4);
  legend1->SetTextAlign(32);
  legend1->SetTextColor(1);
  legend1->SetTextSize(0.04);
  legend1->AddEntry("hi_stasim_eff_pt",l2string,"pl");
  legend1->AddEntry("hi_glbsim_eff_pt",l3string,"pl");
  legend1 ->Draw();
  c1->Write();
  // 
  cout << "Debug 591" << endl;
  TCanvas* c1a = new TCanvas("algEff","Algo Efficiency pt",10,10,700,500);
  c1a->SetFillColor(0);
  c1a->SetGrid(1);
  c1a->SetTicky();
  c1a->SetRightMargin(0.03);
  c1a->SetTopMargin(0.02);
  c1a->cd(); 
  
  hi_glbsta_eff_pt->SetXTitle("p_{T}^{#mu}");
  hi_glbsta_eff_pt->SetYTitle("Efficiency");
  hi_glbsta_eff_pt->SetTitleOffset(1.1,"x");
  hi_glbsta_eff_pt->SetStats(false);

  //hi_glbsta_eff_pt->SetMaximum(1.01);
  //hi_glbsta_eff_pt->SetMinimum(0.8);
  hi_glbsta_eff_pt->SetLineColor(2);
  hi_glbsta_eff_pt->SetMarkerStyle(22);
  hi_glbsta_eff_pt ->DrawCopy("E");
  c1a->Write();
  cout << "Debug 611" << endl;

  TCanvas* c2 = new TCanvas("eff_eta","Efficiency eta",10,10,700,500);
  c2->SetFillColor(0);
  c2->SetGrid(1);
  c2->SetTicky();
  c2->SetRightMargin(0.03);
  c2->SetTopMargin(0.02);
  c2->cd(); 
  cout << "Debug 620" << endl;  
  hi_stasim_eff_eta->SetXTitle("#eta^{#mu}");
  hi_stasim_eff_eta->SetYTitle("Efficiency");
  hi_stasim_eff_eta->SetTitleOffset(1.1,"x");
  hi_stasim_eff_eta->SetTitleOffset(1.15,"y");
  //hi_stasim_eff_eta->SetMaximum(1.02);
  //hi_stasim_eff_eta->SetMinimum(0.5);
  hi_stasim_eff_eta->SetStats(false);
  cout << "Debug 628" << endl;
  hi_stasim_eff_eta->SetLineWidth(1);
  hi_glbsim_eff_eta->SetLineWidth(1);
  hi_stasim_eff_eta->SetLineColor(2);
  hi_glbsim_eff_eta->SetLineColor(4);
  hi_stasim_eff_eta->SetLineStyle(1);
  hi_glbsim_eff_eta->SetLineStyle(1);

  //hi_stasim_eff_eta->SetFillColor(5);
  //hi_glbsim_eff_eta->SetFillColor(7);

  hi_stasim_eff_eta->SetMarkerStyle(21);
  hi_glbsim_eff_eta->SetMarkerStyle(26);
  hi_stasim_eff_eta->SetMarkerColor(2);
  hi_glbsim_eff_eta->SetMarkerColor(4);
  cout << "Debug 643" << endl;
  hi_stasim_eff_eta ->DrawCopy("PE");
  hi_glbsim_eff_eta ->DrawCopy("PEsame");
  //hi_stasim_eff_eta ->DrawCopy("AxisSame");
  cout << "Debug 647" << endl;
  TLegend* legend2 = new TLegend(0.6,0.2,0.8,0.4);
  legend2->SetTextAlign(32);
  legend2->SetTextColor(1);
  legend2->SetTextSize(0.04);
  legend2->AddEntry("hi_stasim_eff_eta",l2string,"pl");
  legend2->AddEntry("hi_glbsim_eff_eta",l3string,"pl");
  legend2 ->Draw();
  c2->Write();
  // 
  cout << "Debug 657" << endl;
  TCanvas* c2a = new TCanvas("algEff_eta","Algo Efficiency eta",10,10,700,500);
  c2a->SetFillColor(0);
  c2a->SetGrid(1);
  c2a->SetTicky();
  c2a->SetRightMargin(0.03);
  c2a->SetTopMargin(0.02);
  c2a->cd(); 
  
  hi_glbsta_eff_eta->SetXTitle("#eta^{#mu}");
  hi_glbsta_eff_eta->SetYTitle("Efficiency");
  hi_glbsta_eff_eta->SetTitleOffset(1.1,"x");
  hi_glbsta_eff_eta->SetStats(false);

  //hi_glbsta_eff_eta->SetMaximum(1.01);
  //hi_glbsta_eff_eta->SetMinimum(0.8);
  hi_glbsta_eff_eta->SetLineColor(2);
  hi_glbsta_eff_eta->SetMarkerStyle(22);
  hi_glbsta_eff_eta ->DrawCopy("E");
  c2a->Write();
  cout << "Debug 677" << endl;

  TCanvas* c3 = new TCanvas("pur","Purity pt",10,10,700,500);
  c3->SetFillColor(0);
  c3->SetGrid(1);
  c3->SetTicky();
  c3->SetRightMargin(0.03);
  c3->SetTopMargin(0.02);
  c3->cd(); 
  cout << "Debug 686" << endl;  
  hi_stasim_pur_pt->SetXTitle("p_{T}^{#mu}");
  hi_stasim_pur_pt->SetYTitle("Purity");
  hi_stasim_pur_pt->SetTitleOffset(1.1,"x");
  hi_stasim_pur_pt->SetTitleOffset(1.15,"y");
  //hi_stasim_pur_pt->SetMaximum(1.02);
  //hi_stasim_pur_pt->SetMinimum(0.5);
  hi_stasim_pur_pt->SetStats(false);
  cout << "Debug 694" << endl;
  hi_stasim_pur_pt->SetLineWidth(1);
  hi_glbsim_pur_pt->SetLineWidth(1);
  hi_stasim_pur_pt->SetLineColor(2);
  hi_glbsim_pur_pt->SetLineColor(4);
  hi_stasim_pur_pt->SetLineStyle(1);
  hi_glbsim_pur_pt->SetLineStyle(1);

  //hi_stasim_pur_pt->SetFillColor(5);
  //hi_glbsim_pur_pt->SetFillColor(7);

  hi_stasim_pur_pt->SetMarkerStyle(21);
  hi_glbsim_pur_pt->SetMarkerStyle(26);
  hi_stasim_pur_pt->SetMarkerColor(2);
  hi_glbsim_pur_pt->SetMarkerColor(4);
  cout << "Debug 709" << endl;
  hi_stasim_pur_pt ->DrawCopy("PE");
  hi_glbsim_pur_pt ->DrawCopy("PEsame");
  //hi_stasim_pur_pt ->DrawCopy("AxisSame");
  cout << "Debug 554" << endl;
  TLegend* legend3 = new TLegend(0.6,0.2,0.8,0.4);
  legend3->SetTextAlign(32);
  legend3->SetTextColor(1);
  legend3->SetTextSize(0.04);
  legend3->AddEntry("hi_stasim_pur_pt",l2string,"pl");
  legend3->AddEntry("hi_glbsim_pur_pt",l3string,"pl");
  legend3 ->Draw();
  c3->Write();
  // 
  cout << "Debug 723" << endl;
  TCanvas* c3a = new TCanvas("algPur","Algo Purity pt",10,10,700,500);
  c3a->SetFillColor(0);
  c3a->SetGrid(1);
  c3a->SetTicky();
  c3a->SetRightMargin(0.03);
  c3a->SetTopMargin(0.02);
  c3a->cd(); 
  
  hi_glbsta_pur_pt->SetXTitle("p_{T}^{#mu}");
  hi_glbsta_pur_pt->SetYTitle("Purity");
  hi_glbsta_pur_pt->SetTitleOffset(1.1,"x");
  hi_glbsta_pur_pt->SetStats(false);

  //hi_glbsta_pur_pt->SetMaximum(1.01);
  //hi_glbsta_pur_pt->SetMinimum(0.8);
  hi_glbsta_pur_pt->SetLineColor(2);
  hi_glbsta_pur_pt->SetMarkerStyle(22);
  hi_glbsta_pur_pt ->DrawCopy("E");
  c3a->Write();
  cout << "Debug 743" << endl;

  TCanvas* c4 = new TCanvas("pur_eta","Purity eta",10,10,700,500);
  c4->SetFillColor(0);
  c4->SetGrid(1);
  c4->SetTicky();
  c4->SetRightMargin(0.03);
  c4->SetTopMargin(0.02);
  c4->cd(); 
  cout << "Debug 752" << endl;  
  hi_stasim_pur_eta->SetXTitle("#eta^{#mu}");
  hi_stasim_pur_eta->SetYTitle("Purity");
  hi_stasim_pur_eta->SetTitleOffset(1.1,"x");
  hi_stasim_pur_eta->SetTitleOffset(1.15,"y");
  //hi_stasim_pur_eta->SetMaximum(1.02);
  //hi_stasim_pur_eta->SetMinimum(0.5);
  hi_stasim_pur_eta->SetStats(false);
  cout << "Debug 760" << endl;
  hi_stasim_pur_eta->SetLineWidth(1);
  hi_glbsim_pur_eta->SetLineWidth(1);
  hi_stasim_pur_eta->SetLineColor(2);
  hi_glbsim_pur_eta->SetLineColor(4);
  hi_stasim_pur_eta->SetLineStyle(1);
  hi_glbsim_pur_eta->SetLineStyle(1);

  //hi_stasim_pur_eta->SetFillColor(5);
  //hi_glbsim_pur_eta->SetFillColor(7);

  hi_stasim_pur_eta->SetMarkerStyle(21);
  hi_glbsim_pur_eta->SetMarkerStyle(26);
  hi_stasim_pur_eta->SetMarkerColor(2);
  hi_glbsim_pur_eta->SetMarkerColor(4);
  cout << "Debug 775" << endl;
  hi_stasim_pur_eta ->DrawCopy("PE");
  hi_glbsim_pur_eta ->DrawCopy("PEsame");
  //hi_stasim_pur_eta ->DrawCopy("AxisSame");
  cout << "Debug 779" << endl;
  TLegend* legend4 = new TLegend(0.6,0.2,0.8,0.4);
  legend4->SetTextAlign(32);
  legend4->SetTextColor(1);
  legend4->SetTextSize(0.04);
  legend4->AddEntry("hi_stasim_pur_eta",l2string,"pl");
  legend4->AddEntry("hi_glbsim_pur_eta",l3string,"pl");
  legend4 ->Draw();
  c4->Write();
  // 
  cout << "Debug 789" << endl;
  TCanvas* c4a = new TCanvas("algPur_eta","Algo Purity eta",10,10,700,500);
  c4a->SetFillColor(0);
  c4a->SetGrid(1);
  c4a->SetTicky();
  c4a->SetRightMargin(0.03);
  c4a->SetTopMargin(0.02);
  c4a->cd(); 
  
  hi_glbsta_pur_eta->SetXTitle("#eta^{#mu}");
  hi_glbsta_pur_eta->SetYTitle("Purity");
  hi_glbsta_pur_eta->SetTitleOffset(1.1,"x");
  hi_glbsta_pur_eta->SetStats(false);

  //hi_glbsta_pur_eta->SetMaximum(1.01);
  //hi_glbsta_pur_eta->SetMinimum(0.8);
  hi_glbsta_pur_eta->SetLineColor(2);
  hi_glbsta_pur_eta->SetMarkerStyle(22);
  hi_glbsta_pur_eta ->DrawCopy("E");
  c4a->Write();
  cout << "Debug 809" << endl;

  hFile->Write();
}

float 
RecoMuonEffAnlyzer::calculateDistance(const math::XYZVector& vect1, const math::XYZVector& vect2) {
  float dEta = vect1.eta() - vect2.eta();
  float dPhi = vect1.phi() - vect2.phi();
  float dPt = sqrt(vect1.perp2()) - sqrt(vect2.perp2());
  float distance = sqrt(pow(dEta,2) + pow(dPhi,2) + pow(dPt,2) );

  return distance;
}

//
// return h1/h2 with recalculated errors
//
TH1F* RecoMuonEffAnlyzer::divideErr(TH1F* h1, TH1F* h2, TH1F* hout) {

  //TH1F* hout = new TH1F(*h1);
  hout->Reset();
  //hout->SetName("DivideErr");
  hout->Divide(h1,h2,1.,1.,"B");

  for (int i = 0; i <= hout->GetNbinsX()+1; i++ ) {
    Float_t tot   = h2->GetBinContent(i) ;
    Float_t tot_e = h2->GetBinError(i);
    Float_t eff = hout->GetBinContent(i) ;
    Float_t Err = 0.;
    if (tot > 0) Err = tot_e / tot * sqrt( eff* (1-eff) );
    if (eff == 1. || isnan(Err) || !isfinite(Err) ) Err=1.e-3;
    hout->SetBinError(i, Err);
  }
  return hout;
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoMuonEffAnlyzer)
