// -*- C++ -*-
//
// Package:    StateOnSurfaceAnalyzer
// Class:      StateOnSurfaceAnalyzer
// 
/**\class StateOnSurfaceAnalyzer StateOnSurfaceAnalyzer.cc UserCode/StateOnSurfaceAnalyzer/src/StateOnSurfaceAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam Everett
//         Created:  Thu Feb  8 18:22:43 CET 2007
// $Id$
//
//


// system include files
#include <memory>
//#include <TFile.h>
#include <TH1.h>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/DiagMatrix.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"



#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/GeomPropagators/interface/StateOnTrackerBound.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

//
// class decleration
//

class StateOnSurfaceAnalyzer : public edm::EDAnalyzer {
public:
  explicit StateOnSurfaceAnalyzer(const edm::ParameterSet&);
  ~StateOnSurfaceAnalyzer();
  
protected:

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual FreeTrajectoryState getFromCLHEP(const Hep3Vector& p3, 
					   const Hep3Vector& r3, 
					   int charge, const HepSymMatrix& cov,
					   const MagneticField* field);
  virtual TrajectoryStateOnSurface propagateSIM(const edm::SimTrackContainer::const_iterator&, const edm::Handle < edm::SimVertexContainer > &);
  virtual double distance(const TrajectoryStateOnSurface&, const TrajectoryStateOnSurface&);  
  // ----------member data ---------------------------
  bool noErrPropMode_;
  //  int trkIndOffset_;
  MuonServiceProxy * theService;
  TrackAssociatorBase * associatorByChi2;
  
  //std::string outputFileName_;
  edm::InputTag simsrc_;
  edm::InputTag tpsrc_;
  edm::InputTag tksrc_;
  edm::InputTag stasrc_;

  TH1F * h_pt, * h_eta, * h_phi;
  TH1F *h_dist1, *h_dist2 , *h_dist3, *h_dist4, *h_dist5, *h_dist6;
  double ptMax_;
};

using namespace edm;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
StateOnSurfaceAnalyzer::StateOnSurfaceAnalyzer(const edm::ParameterSet& iConfig) :
  simsrc_( iConfig.getUntrackedParameter<edm::InputTag>( "simsrc" ) ),
  tpsrc_( iConfig.getUntrackedParameter<edm::InputTag>( "tpsrc" ) ),
  tksrc_( iConfig.getUntrackedParameter<edm::InputTag>( "tksrc" ) ),
  stasrc_( iConfig.getUntrackedParameter<edm::InputTag>( "stasrc" ) ),
  ptMax_( iConfig.getUntrackedParameter<double>( "ptMax" ) )
{
   //now do what ever initialization is needed
 theService = new MuonServiceProxy( iConfig.getParameter<edm::ParameterSet>("ServiceParameters"));

 //  trkIndOffset_ = iConfig.getParameter<int>("trkIndOffset");
  noErrPropMode_ = iConfig.getParameter<bool>("noErrorPropagationMode");

}


StateOnSurfaceAnalyzer::~StateOnSurfaceAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
StateOnSurfaceAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  
  LogDebug("Analyzer") << "++++++++++ Next Event ++++++++++";
  
  theService->update(iSetup);

  /*
   Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel( src_, tracks );
   for( reco::TrackCollection::const_iterator t = tracks->begin(); t != tracks->end(); ++ t ) {
     double pt = t->pt(), eta = t->eta(), phi = t->phi();
     h_pt->Fill( pt );
     h_eta->Fill( eta );
     h_phi->Fill( phi );
   }
   */

  Handle<reco::TrackCollection> stamuons;
  iEvent.getByLabel( stasrc_, stamuons );
  const reco::TrackCollection staTC = *(stamuons.product() );

  Handle<reco::TrackCollection> tkmuons;
  iEvent.getByLabel( tksrc_, tkmuons );
  const reco::TrackCollection tkTC = *(tkmuons.product());

  Handle<SimTrackContainer> simmuons;
  iEvent.getByLabel( simsrc_, simmuons );
  const SimTrackContainer simTC = *(simmuons.product());

  Handle<SimVertexContainer> simVertices;
  iEvent.getByType<SimVertexContainer>(simVertices);

  Handle<TrackingParticleCollection> trackingparticles;
  iEvent.getByLabel( tpsrc_, trackingparticles );

  //make association
  reco::SimToRecoCollection assocTk =
    associatorByChi2->associateSimToReco(tkmuons,trackingparticles,&iEvent);

  reco::SimToRecoCollection assocSTA =
    associatorByChi2->associateSimToReco(stamuons,trackingparticles,&iEvent);

  TrajectoryStateTransform tsTransform;  
  StateOnTrackerBound fromInside(&*theService->propagator("SmartPropagator"));
  StateOnTrackerBound fromOutside(&*theService->propagator("SmartPropagatorOpposite"));

  SimTrackContainer::size_type i=0;
  for(SimTrackContainer::const_iterator simTrack=simTC.begin(); simTrack != simTC.end(); ++simTrack, ++i) {
    TrackingParticleRef tp (trackingparticles,i);
    if( fabs(tp->pdgId()) != 13 ) continue;

    //propagate SIM
    TrajectoryStateOnSurface simTsos = propagateSIM(simTrack,simVertices);
    
    //check Tk
    //    try {
    //      std::vector<std::pair<TrackRef, double> > trackV = assocTk[tp];
    //      for(std::vector<std::pair<TrackRef, double> >::const_iterator it=trackV.begin(); it != trackV.end(); ++it ) {
    for( reco::TrackCollection::const_iterator it=tkTC.begin(); it != tkTC.end(); ++it ) {
	
      //ADAM *it or it->first
      
      TrajectoryStateOnSurface outerTkTsos = tsTransform.outerStateOnSurface(*it,*theService->trackingGeometry(),&*theService->magneticField());
      
	TransientTrack tkTT(*it,&*theService->magneticField(),theService->trackingGeometry());
	TrajectoryStateOnSurface initialTkFts = tkTT.impactPointState();
	
	//propagate Tk
	TrajectoryStateOnSurface tkTsos1 = fromInside(initialTkFts);
	TrajectoryStateOnSurface tkTsos2 = fromInside(outerTkTsos);
	
	//calculae distance
	double distance1 = distance(simTsos,tkTsos1);
	h_dist1->Fill(distance1);
	double distance2 = distance(simTsos,tkTsos2);
	h_dist2->Fill(distance2);
      }
      //    } catch (Exception event) {
      //      LogDebug("Analyzer") << " No Match Found" ;
      //    }
    
    //check STA
      //    try {
      //      std::vector<std::pair<TrackRef, double> > staV = assocSTA[tp];
      //      for(std::vector<std::pair<TrackRef, double> >::const_iterator it=staV.begin(); it != staV.end(); ++it ) {
      for(reco::TrackCollection::const_iterator it=staTC.begin(); it != staTC.end(); ++it ) {
	
	TransientTrack muTT(*it,&*theService->magneticField(),theService->trackingGeometry());
	TrajectoryStateOnSurface ipsMuTsos = muTT.impactPointState();

	TrajectoryStateOnSurface innerMuTsos = tsTransform.innerStateOnSurface(*it,*theService->trackingGeometry(),&*theService->magneticField());

	//propagate STA
	TrajectoryStateOnSurface muTsos1 = fromInside(ipsMuTsos);
	TrajectoryStateOnSurface muTsos2 = fromOutside(innerMuTsos);

	//calculae distance
	double distance3 = distance(simTsos,muTsos1);
	h_dist3->Fill(distance3);
	double distance4 = distance(simTsos,muTsos2);
	h_dist4->Fill(distance4);
      }
      //    } catch (Exception event) {
      //      LogDebug("Analyzer") << " No Match Found" ;
      //    }

  }

  //check Tk - STA
  for(reco::TrackCollection::const_iterator it=staTC.begin(); it != staTC.end(); ++it ) {
    
    TransientTrack muTT(*it,&*theService->magneticField(),theService->trackingGeometry());
    TrajectoryStateOnSurface ipsMuTsos = muTT.impactPointState();
    
    TrajectoryStateOnSurface innerMuTsos = tsTransform.innerStateOnSurface(*it,*theService->trackingGeometry(),&*theService->magneticField());
    
    //propagate STA
    TrajectoryStateOnSurface muTsos1 = fromInside(ipsMuTsos);
    TrajectoryStateOnSurface muTsos2 = fromOutside(innerMuTsos);
    for( reco::TrackCollection::const_iterator it=tkTC.begin(); it != tkTC.end(); ++it ) {
      
      TrajectoryStateOnSurface outerTkTsos = tsTransform.outerStateOnSurface(*it,*theService->trackingGeometry(),&*theService->magneticField());
      
      TransientTrack tkTT(*it,&*theService->magneticField(),theService->trackingGeometry());
      TrajectoryStateOnSurface initialTkFts = tkTT.impactPointState();
      
      //propagate Tk
      TrajectoryStateOnSurface tkTsos1 = fromInside(initialTkFts);
      TrajectoryStateOnSurface tkTsos2 = fromInside(outerTkTsos);
      
      double distance5 = distance(muTsos1,tkTsos2);   
      double distance6 = distance(muTsos2,tkTsos2);   
      
      h_dist5->Fill(distance5);
      h_dist6->Fill(distance6);
    }
  }


}


TrajectoryStateOnSurface StateOnSurfaceAnalyzer::propagateSIM(const SimTrackContainer::const_iterator& tracksCI, const Handle<SimVertexContainer>& simVertices)
{

  int trkPDG = tracksCI->type();
  if (abs(trkPDG) != 13 ) {      
    LogDebug("Analyzer")<<"Skip "<<trkPDG<<std::endl;      
    return TrajectoryStateOnSurface();
  }
  Hep3Vector p3T = tracksCI->momentum().vect();
  if (p3T.mag()< 2.) return TrajectoryStateOnSurface();
  
  int vtxInd = tracksCI->vertIndex();
  //uint trkInd = tracksCI->genpartIndex() - trkIndOffset_;
  Hep3Vector r3T(0.,0.,0.);
  if (vtxInd < 0){
    LogDebug("Analyzer")<<"Track with no vertex, defaulting to (0,0,0)"<<std::endl;
  } else {
    r3T = (*simVertices)[vtxInd].position().vect()*0.1; 
    //seems to be stored in mm --> convert to cm
  }
  HepSymMatrix covT = noErrPropMode_ ? HepSymMatrix(1,1) : HepSymMatrix(6,1); 
  covT *= 1e-20; // initialize to sigma=1e-10 .. should get overwhelmed by MULS
  
  Hep3Vector p3F,r3F; //propagated state
  HepSymMatrix covF(6,0);
  int charge = trkPDG > 0 ? -1 : 1; //works for muons
  
  FreeTrajectoryState ftsTrack = getFromCLHEP(p3T, r3T, charge, covT, &*theService->magneticField());
  FreeTrajectoryState ftsStart = ftsTrack;
  
  StateOnTrackerBound tracker(&*theService->propagator("SmartPropagator"));
  TrajectoryStateOnSurface simTsosFromSimFTS  = tracker(ftsTrack);
  
  return simTsosFromSimFTS;
    
}


// ------------ method called once each job just before starting event loop  ------------
void 
StateOnSurfaceAnalyzer::beginJob(const edm::EventSetup& iSetup)
{
  edm::Service<TFileService> fs;
  h_pt    = fs->make<TH1F>( "pt"  , "p_{t}", 100,  0., ptMax_ );
  h_eta   = fs->make<TH1F>( "eta" , "#eta" , 50, -3, 3 );
  h_phi   = fs->make<TH1F>( "phi" , "#phi" , 50, -M_PI, M_PI );

  h_dist1 = fs->make<TH1F>("h_dist1","Distance SIM Tk_{initial}",100,0,0.1);
  h_dist2 = fs->make<TH1F>("h_dist2","Distance SIM Tk_{outer}",100,0,0.1);
  h_dist3 = fs->make<TH1F>("h_dist3","Distance SIM STA_{initial}",200,0,1);
  h_dist4 = fs->make<TH1F>("h_dist4","Distance SIM STA_{inner}",300,0,1);
  h_dist5 = fs->make<TH1F>("h_dist5","Distance Tk_{outer} STA_{initial}",300,0,1);
  h_dist6 = fs->make<TH1F>("h_dist6","Distance Tk_{outer} STA_{inner}",300,0,1);

  edm::ESHandle<TrackAssociatorBase> theChiAssociator;
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByChi2",theChiAssociator);
  associatorByChi2 = (TrackAssociatorBase *) theChiAssociator.product();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
StateOnSurfaceAnalyzer::endJob() {
}


FreeTrajectoryState
StateOnSurfaceAnalyzer::getFromCLHEP(const Hep3Vector& p3, 
				     const Hep3Vector& r3,
				     int charge, const HepSymMatrix& cov,
				     const MagneticField* field)
{
  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);
  
  CartesianTrajectoryError tCov(cov);
  
  return cov.num_row() == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

double StateOnSurfaceAnalyzer::distance(const TrajectoryStateOnSurface& Tsos1, 
					const TrajectoryStateOnSurface& Tsos2) {
  if( !Tsos1.isValid() || !Tsos2.isValid()) return -1.;
  
  Hep3Vector vect1(Tsos1.globalPosition().x(),Tsos1.globalPosition().y(),Tsos1.globalPosition().z());
  Hep3Vector vect2(Tsos2.globalPosition().x(),Tsos2.globalPosition().y(),Tsos2.globalPosition().z());
  
  return vect1.deltaR(vect2);
}

//define this as a plug-in
DEFINE_FWK_MODULE(StateOnSurfaceAnalyzer);
