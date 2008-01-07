// -*- C++ -*-
// 
// Package:    MatchStudy
// Class:      MatchStudy
// 
/**\class MatchStudy MatchStudy.cc UserCode/MatchStudy/src/MatchStudy.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  A. Everett - Purdue University
//         Created:  Tue Oct  2 12:38:18 EDT 2007
// $Id: MatchStudy.cc,v 1.5 2007/10/20 17:51:45 aeverett Exp $
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/InputTag.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/GeomPropagators/interface/StateOnTrackerBound.h"
#include "TrackingTools/GeomPropagators/interface/StateOnMuonBound.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/GeometrySurface/interface/TangentPlane.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>

using namespace edm;
using namespace std;
using namespace reco;

//
// class decleration
//

class MatchStudy : public edm::EDAnalyzer {
public:
  explicit MatchStudy(const edm::ParameterSet&);
  ~MatchStudy();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  double match_R_IP(const reco::Track&, const reco::Track&) const;
  double match_D(const TrajectoryStateOnSurface&, const TrajectoryStateOnSurface&) const;
  double match_Rmom(const TrajectoryStateOnSurface&, const TrajectoryStateOnSurface&) const;
  double match_Rpos(const TrajectoryStateOnSurface&, const TrajectoryStateOnSurface&) const;
  std::pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface> convertToTSOSTk(const reco::Track&,const reco::Track& ) const;
  std::pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface> convertToTSOSMu(const reco::Track&,const reco::Track& ) const;
  std::pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface> convertToTSOSMuHit(const reco::Track&,const reco::Track& ) const;
  bool samePlane(const TrajectoryStateOnSurface&,const TrajectoryStateOnSurface&) const;
  double matchChiAtSurface(const TrajectoryStateOnSurface& , const TrajectoryStateOnSurface& ) const;
  void FillMtch(const reco::Track&,const reco::Track&,int) const;

      // ----------member data ---------------------------
  TH1S *hEn;
  TH1F *Sim_eta,*Sim_P,*Sim_Pt;
  TH1F **r_ip;
  TH1F **r_tk_mom, **r_tk_pos;
  TH1F **r_mu_mom, **r_mu_pos;
  TH1F **r_muHit_mom, **r_muHit_pos;
  TH1F **d_mu, **d_tk, **d_muHit;
  TH1F **chi2_tk, **surface_tk;
  TH1F **chi2_mu, **surface_mu;
  TH1F **chi2_muHit, **surface_muHit;
  TH1F **matchMethod;

  double theMinP, theMinPt, theMaxChi2, theDeltaEta, theDeltaPhi, theDeltaD, theDeltaR;
  
  const TrackAssociatorBase *tkAssociator_, *muAssociator_;
  std::string tkAssociatorName_, muAssociatorName_;
  edm::InputTag tkName_, tpName_, glbName_, staName_;

  std::string theInPropagatorName;
  std::string theOutPropagatorName;
  MuonServiceProxy *theService;

};

//
// constants, enums and typedefs
//


enum matches {good,bad,wrong};
const Int_t nbins=sizeof(matches)-1;



//
// static data member definitions
//

//
// constructors and destructor
//
MatchStudy::MatchStudy(const edm::ParameterSet& iConfig)

{
  // service parameters
  ParameterSet serviceParameters = iConfig.getParameter<ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  
  //now do what ever initialization is needed
  tkAssociatorName_ = iConfig.getUntrackedParameter<std::string>("tkAssociator");
  muAssociatorName_ = iConfig.getUntrackedParameter<std::string>("muAssociator");

  tpName_ = iConfig.getUntrackedParameter<edm::InputTag>("tpLabel");
  tkName_ = iConfig.getUntrackedParameter<edm::InputTag>("tkLabel");
  staName_ = iConfig.getUntrackedParameter<edm::InputTag>("muLabel");
  glbName_ = iConfig.getUntrackedParameter<edm::InputTag>("glbLabel");


  theOutPropagatorName = iConfig.getParameter<string>("StateOnTrackerBoundOutPropagator");
  theInPropagatorName = iConfig.getParameter<string>("StateOnMuonBoundInPropagator");

  theMaxChi2  = iConfig.getParameter<double>("Chi2Cut");
  theDeltaEta = iConfig.getParameter<double>("DeltaEtaCut");
  theDeltaPhi = iConfig.getParameter<double>("DeltaPhiCut");
  theDeltaD   = iConfig.getParameter<double>("DeltaDCut");
  theDeltaR   = iConfig.getParameter<double>("DeltaRCut");
  theMinP     = iConfig.getParameter<double>("MinP");
  theMinPt    = iConfig.getParameter<double>("MinPt");
}


MatchStudy::~MatchStudy()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete theService;
}


//
// member functions
//


// ------------ method called once each job just before starting event loop  ------------



void 
MatchStudy::beginJob(const edm::EventSetup& setup)
{
  // Tk Associator
  edm::ESHandle<TrackAssociatorBase> tkassociatorHandle;
  setup.get<TrackAssociatorRecord>().get(tkAssociatorName_,tkassociatorHandle);
  tkAssociator_ = tkassociatorHandle.product();

  // Mu Associator
  edm::ESHandle<TrackAssociatorBase> muassociatorHandle;
  setup.get<TrackAssociatorRecord>().get(muAssociatorName_,muassociatorHandle);
  muAssociator_ = muassociatorHandle.product();


  edm::Service<TFileService> fs;


  hEn = fs->make<TH1S>("hEn","iEvent.id().event() Entries",1000,0,1000);
  Sim_eta = fs->make<TH1F>("Sim_eta","simulated eta distribution",100,0,2.5);
  Sim_P = fs->make<TH1F>("Sim_P","simulated P distribution",5000,0,5000);
  Sim_Pt = fs->make<TH1F>("Sim_Pt","simulated Pt distribution",5000,0,5000);

  r_ip = new TH1F*[nbins];
  r_tk_mom = new TH1F*[nbins];
  r_tk_pos = new TH1F*[nbins];
  r_mu_mom = new TH1F*[nbins];
  r_mu_pos = new TH1F*[nbins];
  r_muHit_mom = new TH1F*[nbins];
  r_muHit_pos = new TH1F*[nbins];
  d_mu = new TH1F*[nbins];
  d_muHit = new TH1F*[nbins];
  d_tk = new TH1F*[nbins];
  chi2_tk = new TH1F*[nbins];
  chi2_mu = new TH1F*[nbins];
  chi2_muHit = new TH1F*[nbins];
  surface_tk = new TH1F*[nbins];
  surface_mu = new TH1F*[nbins];
  surface_muHit = new TH1F*[nbins];
  matchMethod = new TH1F*[nbins];


  // label histograms according to enumeration 
  for (Int_t i=0;i<nbins;i++) {
    
    std::stringstream c;
    c << i;
    TString en(c.str());

    r_ip[i] = fs->make<TH1F>(TString("r_ip")+en,"matched R at IP",100,0.,1.);
    r_tk_mom[i] = fs->make<TH1F>(TString("r_tk_mom")+en,"R_{mom} at Tk Surface",100,0.,1.);
    r_tk_pos[i] = fs->make<TH1F>(TString("r_tk_pos")+en,"R_{pos} at Tk Surface",100,0.,1.);
    r_mu_mom[i] = fs->make<TH1F>(TString("r_mu_mom")+en,"R_{mom} at Mu Surface",100,0.,1.);
    r_mu_pos[i] = fs->make<TH1F>(TString("r_mu_pos")+en,"R_{pos} at Mu Surface",100,0.,1.);
    r_muHit_mom[i] = fs->make<TH1F>(TString("r_muHit_mom")+en,"R_{mom} at Mu Hit Surface",100,0.,1.);
    r_muHit_pos[i] = fs->make<TH1F>(TString("r_muHit_pos")+en,"R_{pos} at Mu Hit Surface",100,0.,1.);
    
    d_mu[i] = fs->make<TH1F>(TString("d_mu")+en,"D at Mu Surface",100,0.,100.);
    d_muHit[i] = fs->make<TH1F>(TString("d_muHit")+en,"D at Mu Hit Surface",100,0.,100.);
    d_tk[i] = fs->make<TH1F>(TString("d_tk")+en,"D at Tk Surface",100,0.,100.);
    
    chi2_tk[i] = fs->make<TH1F>(TString("chi2_tk")+en,"#chi^{2} of all tracks on Tk Surface",501,-2.,1000.);
    chi2_mu[i] = fs->make<TH1F>(TString("chi2_mu")+en,"#chi^{2} of all tracks on Mu Surface",501,-2.,1000.);
    chi2_muHit[i] = fs->make<TH1F>(TString("chi2_muHit")+en,"#chi^{2} of all tracks on Mu Hit Surface",501,-2.,1000.);
    
    surface_tk[i] = fs->make<TH1F>(TString("surface_tk")+en,"Pass/Fail Tk Surface",3,-0.5,2.5);
    surface_mu[i] = fs->make<TH1F>(TString("surface_mu")+en,"Pass/Fail Mu Surface",3,-0.5,2.5);
    surface_muHit[i] = fs->make<TH1F>(TString("surface_muHit")+en,"Pass/Fail Mu Hit Surface",3,-0.5,2.5);
    
    matchMethod[i] = fs->make<TH1F>(TString("matchMethod")+en,"MatchMethod",11,-0.5,10.5);
  }
 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MatchStudy::endJob() {
}



// ------------ method called to for each event  ------------
void
MatchStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  using namespace edm;
  using namespace std;
  using namespace reco;

  hEn->Fill(iEvent.id().event());

  theService->update(iSetup);

  Handle<TrackingParticleCollection> tpHandle;
  iEvent.getByLabel(tpName_,tpHandle);
  const TrackingParticleCollection tpColl = *(tpHandle.product());

  Handle<reco::TrackCollection> staHandle;
  iEvent.getByLabel(staName_,staHandle);
  const reco::TrackCollection staColl = *(staHandle.product());

  Handle<reco::TrackCollection> tkHandle;
  iEvent.getByLabel(tkName_,tkHandle);
  const reco::TrackCollection tkColl = *(tkHandle.product());  

  reco::RecoToSimCollection tkrecoToSimCollection;
  reco::SimToRecoCollection tksimToRecoCollection;
  
  tkrecoToSimCollection = tkAssociator_->associateRecoToSim(tkHandle,tpHandle,&iEvent);
  tksimToRecoCollection = tkAssociator_->associateSimToReco(tkHandle,tpHandle,&iEvent);
  
  reco::RecoToSimCollection starecoToSimCollection;
  reco::SimToRecoCollection stasimToRecoCollection;
  starecoToSimCollection = muAssociator_->associateRecoToSim(staHandle,tpHandle,&iEvent);
  stasimToRecoCollection = muAssociator_->associateSimToReco(staHandle,tpHandle,&iEvent);

  //  LogDebug("MatchStudy") << "Association size: tk2sim =" << tkrecoToSimCollection.size() << ", sim2tk =" << tksimToRecoCollection.size();
  //  LogDebug("MatchStudy") << "Association size: sta2sim =" << starecoToSimCollection.size() << ", sim2sta =" << stasimToRecoCollection.size();
  
  for (TrackingParticleCollection::size_type i=0; i<tpColl.size(); ++i){
    TrackingParticleRef tp(tpHandle,i);
    
    // Get Momenta and Eta of simulated muons 
    if(abs(tp->pdgId())==13){
      Sim_Pt->Fill(tp->pt());
      Sim_P->Fill(tp->p());
      Sim_eta->Fill(tp->eta());
    }
    
    std::vector<std::pair<reco::TrackRef, double> > rvSta;
    reco::TrackRef iSta;
    if(stasimToRecoCollection.find(tp) != stasimToRecoCollection.end()){
      rvSta = stasimToRecoCollection[tp];
      if(rvSta.size() != 0) {
	iSta = rvSta.begin()->first;
      }
    }
    
    std::vector<std::pair<reco::TrackRef, double> > rvTk;
    reco::TrackRef iTk;
    if(tksimToRecoCollection.find(tp) != tksimToRecoCollection.end()){
      rvTk = tksimToRecoCollection[tp];
      if(rvTk.size() != 0) {
	iTk = rvTk.begin()->first;
      }
    }
    
    // match to correct muon in tracker  
    FillMtch(*iSta,*iTk,good); 
    
    for( reco::TrackCollection::size_type it=0; it<tkColl.size(); it++){

      reco::TrackRef T(tkHandle,it);
      
      if(T != iTk){ 
	//  match to any tracker track
	FillMtch(*iSta,*T,bad); 
	// match to wrong muon in the tracker  
	if(tkrecoToSimCollection.find(T)!= tkrecoToSimCollection.end())
	  FillMtch(*iSta,*T,wrong);
      }
    }
    
  }
}


// ------------------- Member Function implementations --------------------------------------

void MatchStudy::FillMtch(const reco::Track& staCand,const reco::Track& tkCand,int en) const {

  std::pair<TrajectoryStateOnSurface, TrajectoryStateOnSurface> tsosPairTk 
    = convertToTSOSTk(staCand,tkCand);
  
  r_ip[en]->Fill( match_R_IP(staCand,tkCand));
	      
  
  bool sameSurface = samePlane(tsosPairTk.first,tsosPairTk.second);
  surface_tk[en]->Fill(sameSurface);
  
  if( sameSurface ) {
    r_tk_mom[en]->Fill( match_Rmom(tsosPairTk.first,tsosPairTk.second) );
    r_tk_pos[en]->Fill( match_Rpos(tsosPairTk.first,tsosPairTk.second) );
    d_tk[en]->Fill( match_D(tsosPairTk.first,tsosPairTk.second) );
  }
  
    double chi2 = -1;
    if( sameSurface ) {
      chi2 = matchChiAtSurface(tsosPairTk.first, tsosPairTk.second);
      chi2_tk[en]->Fill(chi2);
    }
    
    ///
    std::pair<TrajectoryStateOnSurface, TrajectoryStateOnSurface> tsosPairMu 
      = convertToTSOSMu(staCand,tkCand);
    
    bool sameSurfaceMu = samePlane(tsosPairMu.first,tsosPairMu.second);
    surface_mu[en]->Fill(sameSurfaceMu);
    
    if( sameSurfaceMu ) {
      r_mu_mom[en]->Fill( match_Rmom(tsosPairMu.first,tsosPairMu.second) );
      r_mu_pos[en]->Fill( match_Rpos(tsosPairMu.first,tsosPairMu.second) );
      d_mu[en]->Fill( match_D(tsosPairMu.first,tsosPairMu.second) );
    }
    
    double chi2Mu = -1;
    if( sameSurfaceMu ) { 
      chi2Mu = matchChiAtSurface(tsosPairMu.first, tsosPairMu.second);
      chi2_mu[en]->Fill(chi2Mu);
    }
    
    ///
    std::pair<TrajectoryStateOnSurface, TrajectoryStateOnSurface> tsosPairMuHit 
      = convertToTSOSMuHit(staCand,tkCand);
    
    bool sameSurfaceMuHit = samePlane(tsosPairMuHit.first,tsosPairMuHit.second);
    surface_muHit[en]->Fill(sameSurfaceMuHit);
    
    if( sameSurfaceMuHit ) {
      r_muHit_mom[en]->Fill( match_Rmom(tsosPairMuHit.first,tsosPairMuHit.second) );
      r_muHit_pos[en]->Fill( match_Rpos(tsosPairMuHit.first,tsosPairMuHit.second) );
      d_muHit[en]->Fill( match_D(tsosPairMuHit.first,tsosPairMuHit.second) );
    }
    
    double chi2MuHit = -1;
    if( sameSurfaceMuHit ) {
      chi2MuHit = matchChiAtSurface(tsosPairMuHit.first, tsosPairMuHit.second);
      chi2_muHit[en]->Fill(chi2MuHit);
    }
    
    /// check match method
    
    matchMethod[en]->Fill(0);
    if(sameSurfaceMuHit) {
      if (match_Rpos(tsosPairMuHit.first,tsosPairMuHit.second) < theDeltaR) {
	matchMethod[en]->Fill(1);
      } else if (match_D(tsosPairMuHit.first,tsosPairMuHit.second) < theDeltaD) {
	matchMethod[en]->Fill(2);
      } else if(matchChiAtSurface(tsosPairMuHit.first,tsosPairMuHit.second) < theMaxChi2) {
	matchMethod[en]->Fill(3); 
      } else {
	matchMethod[en]->Fill(10);
      }
    }
    
}


std::pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>
MatchStudy::convertToTSOSTk(const reco::Track& staCand,
			    const reco::Track& tkCand) const {
  
  const string category = "MatchStudy";
  
  TransientTrack muTT(staCand,&*theService->magneticField(),theService->trackingGeometry());
  TrajectoryStateOnSurface impactMuTSOS = muTT.impactPointState();

  TrajectoryStateOnSurface outerTkTsos;

  // make sure the tracker Track has enough momentum to reach the muon chambers
  if ( !(tkCand.p() < theMinP || tkCand.pt() < theMinPt )) {
    TrajectoryStateTransform tsTransform;
    outerTkTsos = tsTransform.outerStateOnSurface(tkCand,*theService->trackingGeometry(),&*theService->magneticField());
  }
  
  if ( !impactMuTSOS.isValid() || !outerTkTsos.isValid() ) return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(impactMuTSOS,outerTkTsos);
  
  // define StateOnTrackerBound objects  
  StateOnTrackerBound fromInside(&*theService->propagator(theOutPropagatorName));
  
  // extrapolate to outer tracker surface
  TrajectoryStateOnSurface tkTsosFromMu = fromInside(impactMuTSOS);
  TrajectoryStateOnSurface tkTsosFromTk = fromInside(outerTkTsos);

    
  if( !samePlane(tkTsosFromMu,tkTsosFromTk)) {
    bool same1, same2;
    //propagate tk to same surface as muon
    TrajectoryStateOnSurface newTkTsosFromTk, newTkTsosFromMu;
    if( tkTsosFromMu.isValid() ) newTkTsosFromTk = theService->propagator(theOutPropagatorName)->propagate(outerTkTsos,tkTsosFromMu.surface());
    same1 =  samePlane(newTkTsosFromTk,tkTsosFromMu);
    LogDebug(category) << "Propagating to same tracker surface (Mu):" << same1;
    if( !same1 ) {
      if( tkTsosFromTk.isValid() ) newTkTsosFromMu = theService->propagator(theOutPropagatorName)->propagate(impactMuTSOS,tkTsosFromTk.surface());
      same2 =  samePlane(newTkTsosFromMu,tkTsosFromTk);
      LogDebug(category) << "Propagating to same tracker surface (Tk):" << same2;
    }
    if(same1) tkTsosFromTk = newTkTsosFromTk;
    else if(same2) tkTsosFromMu = newTkTsosFromMu;
    else  {
      LogDebug(category) << "Could not propagate Muon and Tracker track to the same tracker bound!";
      return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(impactMuTSOS, outerTkTsos);
    }
  }
  
  return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(tkTsosFromMu, tkTsosFromTk);
}

std::pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>
MatchStudy::convertToTSOSMu(const reco::Track& staCand,
			    const reco::Track& tkCand) const {
  
  const string category = "MatchStudy";
  
  TransientTrack muTT(staCand,&*theService->magneticField(),theService->trackingGeometry());
  //TrajectoryStateOnSurface impactMuTSOS = muTT.impactPointState();
  TrajectoryStateOnSurface innerMuTSOS = muTT.innermostMeasurementState();

  TrajectoryStateOnSurface outerTkTsos;

  // make sure the tracker Track has enough momentum to reach the muon chambers
  if ( !(tkCand.p() < theMinP || tkCand.pt() < theMinPt )) {
    TrajectoryStateTransform tsTransform;
    outerTkTsos = tsTransform.outerStateOnSurface(tkCand,*theService->trackingGeometry(),&*theService->magneticField());
  }
  
  if ( !innerMuTSOS.isValid() || !outerTkTsos.isValid() ) return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(innerMuTSOS,outerTkTsos);
  
  // define StateOnTrackerBound objects  
  StateOnMuonBound fromInside(&*theService->propagator(theOutPropagatorName));
  StateOnMuonBound fromOutside(&*theService->propagator(theInPropagatorName));
  
  // extrapolate to outer tracker surface
  TrajectoryStateOnSurface muTsosFromMu = fromOutside(innerMuTSOS);
  TrajectoryStateOnSurface muTsosFromTk = fromInside(outerTkTsos);

    
  if( !samePlane(muTsosFromMu,muTsosFromTk)) {
    bool same1, same2;
    //propagate tk to same surface as muon
    TrajectoryStateOnSurface newMuTsosFromTk, newMuTsosFromMu;
    if( muTsosFromMu.isValid() ) newMuTsosFromTk = theService->propagator(theOutPropagatorName)->propagate(outerTkTsos,muTsosFromMu.surface());
    same1 =  samePlane(newMuTsosFromTk,muTsosFromMu);
    LogDebug(category) << "Propagating to same muon surface (Mu):" << same1;
    if( !same1 ) {
      if( muTsosFromTk.isValid() ) newMuTsosFromMu = theService->propagator(theInPropagatorName)->propagate(innerMuTSOS,muTsosFromTk.surface());
      same2 =  samePlane(newMuTsosFromMu,muTsosFromTk);
      LogDebug(category) << "Propagating to same muon surface (Tk):" << same2;
    }
    if(same1) muTsosFromTk = newMuTsosFromTk;
    else if(same2) muTsosFromMu = newMuTsosFromMu;
    else {
      LogDebug(category) << "Could not propagate Muon and Tracker track to the same muon bound!";
     //  ReferenceCountingPointer<TangentPlane> p1(newMuTsosFromTk.surface().tangentPlane(innerMuTSOS.localPosition()));
//       ReferenceCountingPointer<TangentPlane> p2(outerMuTSOS.surface().tangentPlane(outerMuTSOS.localPosition()));
//       LogDebug(category) << "tilt = " p1->normalVector().dot(p2->normalVector()) << ", dist =" << p1->toLocal(p2->position())).z() << ;
    return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(innerMuTSOS, outerTkTsos);
  }
}
  
  return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(muTsosFromMu, muTsosFromTk);
}


std::pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>
MatchStudy::convertToTSOSMuHit(const reco::Track& staCand,
			       const reco::Track& tkCand) const {
  
  const string category = "MatchStudy";
  
  TransientTrack muTT(staCand,&*theService->magneticField(),theService->trackingGeometry());
  TrajectoryStateOnSurface innerMuTSOS = muTT.innermostMeasurementState();

  TrajectoryStateOnSurface outerTkTsos;

  // make sure the tracker Track has enough momentum to reach the muon chambers
  if ( !(tkCand.p() < theMinP || tkCand.pt() < theMinPt )) {
    TrajectoryStateTransform tsTransform;
    outerTkTsos = tsTransform.outerStateOnSurface(tkCand,*theService->trackingGeometry(),&*theService->magneticField());
  }
  
  if ( !innerMuTSOS.isValid() || !outerTkTsos.isValid() ) return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(innerMuTSOS,outerTkTsos);
  
  TrajectoryStateOnSurface tkAtMu = theService->propagator(theOutPropagatorName)->propagate(outerTkTsos,theService->trackingGeometry()->idToDet( DetId(staCand.innerDetId()) )->surface());
  //TrajectoryStateOnSurface tkAtMu = theService->propagator(theOutPropagatorName)->propagate(outerTkTsos,innerMuTSOS.surface());
  
  
  if( !samePlane(innerMuTSOS,tkAtMu)) {
    LogDebug(category) << "Could not propagate Muon and Tracker track to the same muon hit surface!";
    return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(innerMuTSOS, outerTkTsos);    
  }
  
  return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(innerMuTSOS, tkAtMu);
}


bool MatchStudy::samePlane(const TrajectoryStateOnSurface& tsos1,
			   const TrajectoryStateOnSurface& tsos2) const
{
  if( !tsos1.isValid() || !tsos2.isValid()) return false;
  const string category = "GlobalMuonTrackMatcher";

  const float maxtilt = 0.999;
  const float maxdist = 0.01; // in cm

  ReferenceCountingPointer<TangentPlane> p1(tsos1.surface().tangentPlane(tsos1.localPosition()));
  ReferenceCountingPointer<TangentPlane> p2(tsos2.surface().tangentPlane(tsos2.localPosition()));

  bool returnValue =  ( (fabs(p1->normalVector().dot(p2->normalVector())) > maxtilt) || (fabs((p1->toLocal(p2->position())).z()) < maxdist) ) ? true : false;

  return returnValue; 
  
}

double 
MatchStudy::matchChiAtSurface(const TrajectoryStateOnSurface& tsos1, 
			      const TrajectoryStateOnSurface& tsos2) const {
  
  const string category = "MatchStudy";
  
  if ( !tsos1.isValid() || !tsos2.isValid() ) return -1.;

  AlgebraicVector5 v(tsos1.localParameters().vector() - tsos2.localParameters().vector());
  AlgebraicSymMatrix55 m(tsos1.localError().matrix() + tsos2.localError().matrix());
  //LogDebug(category) << "vector v " << v;

  int ierr = ! m.Invert();

  if (ierr != 0) edm::LogInfo(category) << "Error inversing covariance matrix";

  double est = ROOT::Math::Similarity(v,m);

  //LogDebug(category) << "Chi2 " << est;

/*
  GlobalVector x = tsos1.globalParameters().position() - tsos2.globalParameters().position();
  AlgebraicVector v1(3); v1[0] = x.x(); v1[1] = x.y(); v1[2] = x.z();
  AlgebraicSymMatrix m1(tsos1.cartesianError().position().matrix() + tsos2.cartesianError().position().matrix());
  m1.invert(ierr);
  double est1 = m1.similarity(v1);
*/

  return est;

}

double
MatchStudy::match_R_IP(const reco::Track& staTrack, const reco::Track& tkTrack) const {
  return (deltaR<double>(staTrack.eta(),staTrack.phi(),
			 tkTrack.eta(),tkTrack.phi()));
}

double
MatchStudy::match_Rmom(const TrajectoryStateOnSurface& sta, const TrajectoryStateOnSurface& tk) const {
  return (deltaR<double>(sta.globalMomentum().eta(),sta.globalMomentum().phi(),
			 tk.globalMomentum().eta(),tk.globalMomentum().phi()));
}

double
MatchStudy::match_Rpos(const TrajectoryStateOnSurface& sta, const TrajectoryStateOnSurface& tk) const {
  return (deltaR<double>(sta.globalPosition().eta(),sta.globalPosition().phi(),
			 tk.globalPosition().eta(),tk.globalPosition().phi()));
}

double
MatchStudy::match_D(const TrajectoryStateOnSurface& sta, const TrajectoryStateOnSurface& tk) const {
  return (sta.globalPosition() - tk.globalPosition()).mag();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MatchStudy);
