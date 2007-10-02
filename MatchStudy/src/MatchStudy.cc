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

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include <TH1.h>
#include <TH2.h>

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
  bool samePlane(const TrajectoryStateOnSurface&,const TrajectoryStateOnSurface&) const;
  double matchChiAtSurface(const TrajectoryStateOnSurface& , const TrajectoryStateOnSurface& ) const;
      // ----------member data ---------------------------
  TH1F *r_ip, *r_tk_mom, *r_tk_pos;
  TH1F        *r_mu_mom, *r_mu_pos;
  TH1F        *d_mu, *d_tk;
  TH1F *chi2_tk_all, *surface_tk;
  TH1F *chi2_mu_all, *surface_mu;

  double theMinP, theMinPt, theMaxChi2, theDeltaEta, theDeltaPhi;
  
  InputTag staTag, tkTag;
  std::string theInPropagatorName;
  std::string theOutPropagatorName;
  MuonServiceProxy *theService;
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
MatchStudy::MatchStudy(const edm::ParameterSet& iConfig)

{
  // service parameters
  ParameterSet serviceParameters = iConfig.getParameter<ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  
   //now do what ever initialization is needed
  staTag = iConfig.getParameter<InputTag>("StaCollection");
  tkTag = iConfig.getParameter<InputTag>("TkCollection");

  theOutPropagatorName = iConfig.getParameter<string>("StateOnTrackerBoundOutPropagator");
  theInPropagatorName = iConfig.getParameter<string>("StateOnMuonBoundInPropagator");

  theMaxChi2 =  iConfig.getParameter<double>("Chi2Cut");
  theDeltaEta = iConfig.getParameter<double>("DeltaEtaCut");
  theDeltaPhi = iConfig.getParameter<double>("DeltaPhiCut");
  theMinP = iConfig.getParameter<double>("MinP");
  theMinPt = iConfig.getParameter<double>("MinPt");
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

// ------------ method called to for each event  ------------
void
MatchStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

  theService->update(iSetup);

  Handle<reco::TrackCollection> staTracksH;
  iEvent.getByLabel(staTag,staTracksH);
  const reco::TrackCollection staTracks = *(staTracksH.product());

  Handle<reco::TrackCollection> tkTracksH;
  iEvent.getByLabel(tkTag,tkTracksH);
  const reco::TrackCollection tkTracks = *(tkTracksH.product());  

  for(reco::TrackCollection::const_iterator iSta=staTracks.begin(); iSta!=staTracks.end();++iSta) {
    for (reco::TrackCollection::const_iterator iTk = tkTracks.begin(); iTk != tkTracks.end(); ++iTk) {
      r_ip->Fill( match_R_IP(*iSta,*iTk) );
      
      std::pair<TrajectoryStateOnSurface, TrajectoryStateOnSurface> tsosPairTk 
	= convertToTSOSTk(*iSta,*iTk);

      r_tk_mom->Fill( match_Rmom(tsosPairTk.first,tsosPairTk.second) );
      r_tk_pos->Fill( match_Rmom(tsosPairTk.first,tsosPairTk.second) );
      d_tk->Fill( match_D(tsosPairTk.first,tsosPairTk.second) );

      bool sameSurface = samePlane(tsosPairTk.first,tsosPairTk.second);
      surface_tk->Fill(sameSurface);

      double chi2 = -1;
      if( sameSurface ) 
	chi2 = matchChiAtSurface(tsosPairTk.first, tsosPairTk.second);
      chi2_tk_all->Fill(chi2);

      ///
      std::pair<TrajectoryStateOnSurface, TrajectoryStateOnSurface> tsosPairMu 
	= convertToTSOSMu(*iSta,*iTk);

      r_mu_mom->Fill( match_Rmom(tsosPairMu.first,tsosPairMu.second) );
      r_mu_pos->Fill( match_Rmom(tsosPairMu.first,tsosPairMu.second) );
      d_mu->Fill( match_D(tsosPairMu.first,tsosPairMu.second) );

      bool sameSurfaceMu = samePlane(tsosPairMu.first,tsosPairMu.second);
      surface_mu->Fill(sameSurfaceMu);

      double chi2Mu = -1;
      if( sameSurfaceMu ) 
	chi2Mu = matchChiAtSurface(tsosPairMu.first, tsosPairMu.second);
      chi2_mu_all->Fill(chi2Mu);

    }
    
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
MatchStudy::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  r_ip = fs->make<TH1F>("r_ip","R at IP",100,0.,1.);

  r_tk_mom = fs->make<TH1F>("r_tk_mom","R_{mom} at Tk Surface",100,0.,1.);
  r_tk_pos = fs->make<TH1F>("r_tk_pos","R_{pos} at Tk Surface",100,0.,1.);

  r_mu_mom = fs->make<TH1F>("r_mu_mom","R_{mom} at Tk Surface",100,0.,1.);
  r_mu_pos = fs->make<TH1F>("r_mu_pos","R_{pos} at Tk Surface",100,0.,1.);

  d_mu = fs->make<TH1F>("d_mu","D at Mu Surface",100,0.,100.);
  d_tk = fs->make<TH1F>("d_tk","D at Tk Surface",100,0.,100.);

  chi2_tk_all = fs->make<TH1F>("chi2_tk_all","#chi^{2} of all tracks on Tk Surface",501,-2.,1000.);
  chi2_mu_all = fs->make<TH1F>("chi2_mu_all","#chi^{2} of all tracks on Mu Surface",501,-2.,1000.);

  surface_tk = fs->make<TH1F>("surface_tk","Pass/Fail Tk Surface",3,-0.5,2.5);
  surface_mu = fs->make<TH1F>("surface_mu","Pass/Fail Mu Surface",3,-0.5,2.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MatchStudy::endJob() {
}

std::pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>
MatchStudy::convertToTSOSTk(const reco::Track& staCand,
			    const reco::Track& tkCand) const {
  
  const string category = "MatchStudy";
  
  LogDebug(category);

  TransientTrack muTT(staCand,&*theService->magneticField(),theService->trackingGeometry());
  TrajectoryStateOnSurface impactMuTSOS = muTT.impactPointState();

  TrajectoryStateOnSurface outerTkTsos;

  LogDebug(category);

  // make sure the tracker Track has enough momentum to reach the muon chambers
  if ( !(tkCand.p() < theMinP || tkCand.pt() < theMinPt )) {
    TrajectoryStateTransform tsTransform;
    outerTkTsos = tsTransform.outerStateOnSurface(tkCand,*theService->trackingGeometry(),&*theService->magneticField());
  }
  
  LogDebug(category);

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
    LogDebug(category) << "Propagating to same surface (Mu):" << same1;
    if( !same1 ) {
      if( tkTsosFromTk.isValid() ) newTkTsosFromMu = theService->propagator(theOutPropagatorName)->propagate(impactMuTSOS,tkTsosFromTk.surface());
      same2 =  samePlane(newTkTsosFromMu,tkTsosFromTk);
      LogDebug(category) << "Propagating to same surface (Tk):" << same2;
    }
    if(same1) tkTsosFromTk = newTkTsosFromTk;
    else if(same2) tkTsosFromMu = newTkTsosFromMu;
    else  LogDebug(category) << "Could not propagate Muon and Tracker track to the same tracker bound!";
  }
  
  return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(tkTsosFromMu, tkTsosFromTk);
}

std::pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>
MatchStudy::convertToTSOSMu(const reco::Track& staCand,
			    const reco::Track& tkCand) const {
  
  const string category = "MatchStudy";
  
  LogDebug(category);

  TransientTrack muTT(staCand,&*theService->magneticField(),theService->trackingGeometry());
  //TrajectoryStateOnSurface impactMuTSOS = muTT.impactPointState();
  TrajectoryStateOnSurface innerMuTSOS = muTT.innermostMeasurementState();

  TrajectoryStateOnSurface outerTkTsos;

  LogDebug(category);

  // make sure the tracker Track has enough momentum to reach the muon chambers
  if ( !(tkCand.p() < theMinP || tkCand.pt() < theMinPt )) {
    TrajectoryStateTransform tsTransform;
    outerTkTsos = tsTransform.outerStateOnSurface(tkCand,*theService->trackingGeometry(),&*theService->magneticField());
  }
  
  LogDebug(category);

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
    LogDebug(category) << "Propagating to same surface (Mu):" << same1;
    if( !same1 ) {
      if( muTsosFromTk.isValid() ) newMuTsosFromMu = theService->propagator(theOutPropagatorName)->propagate(innerMuTSOS,muTsosFromTk.surface());
      same2 =  samePlane(newMuTsosFromMu,muTsosFromTk);
      LogDebug(category) << "Propagating to same surface (Tk):" << same2;
    }
    if(same1) muTsosFromTk = newMuTsosFromTk;
    else if(same2) muTsosFromMu = newMuTsosFromMu;
    else  LogDebug(category) << "Could not propagate Muon and Tracker track to the same tracker bound!";
  }
  
  return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(muTsosFromMu, muTsosFromTk);
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
  LogDebug(category) << "vector v " << v;

  int ierr = ! m.Invert();

  if (ierr != 0) edm::LogInfo(category) << "Error inversing covariance matrix";

  double est = ROOT::Math::Similarity(v,m);

  LogDebug(category) << "Chi2 " << est;

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
