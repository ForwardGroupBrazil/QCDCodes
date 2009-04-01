// -*- C++ -*-
//
// Package:    GlbSelectorStudy
// Class:      GlbSelectorStudy
// 
/**\class GlbSelectorStudy GlbSelectorStudy.cc UserCode/GlbSelectorStudy/src/GlbSelectorStudy.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam A Everett
//         Created:  Tue Mar 31 16:20:19 EDT 2009
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "TrackingTools/GeomPropagators/interface/TrackerBounds.h"
#include "TrackingTools/GeomPropagators/interface/MuonBounds.h"

#include "UserCode/GlbSelectorStudy/interface/IDconverttoBinNum.h"
#include "UserCode/GlbSelectorStudy/interface/MotherSearch.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"

#include "TH1.h"
#include "TFile.h"

//
// class decleration
//

class GlbSelectorStudy : public edm::EDAnalyzer {
   public:
      explicit GlbSelectorStudy(const edm::ParameterSet&);
      ~GlbSelectorStudy();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  bool inTrackerBound_;
  bool inMuonBound_;
  bool inCal_;
  bool simMuon_;

  bool doAssoc_;

  std::string theCategory;

  IDconverttoBinNum wantMotherBin;

  //edm::InputTag simLabel_;
  //edm::InputTag trackingParticleLabel_;
  //edm::InputTag glbMuLabel_;

  //edm::InputTag glbMuAssocLabel_;

  edm::InputTag simLabel_;
  edm::InputTag trkMuLabel_;
  edm::InputTag staMuLabel_;
  edm::InputTag glbMuLabel_;
  edm::InputTag muonLabel_;

  edm::InputTag trkMuAssocLabel_;
  edm::InputTag staMuAssocLabel_;
  edm::InputTag glbMuAssocLabel_;

  TrackAssociatorBase* trkMuAssociator_, * staMuAssociator_, * glbMuAssociator_;

  //TFileDirectory * prompt; //, *promptSta, *promptTk, *promptGlb;
  //TFileDirectory decay;//, *decaySta, *decayTk, *decayGlb;
  //TFileDirectory outside;//, *outsideSta, *outsideTk, *outsideGlb;

  struct MuonME;
  MuonME *promptME_,*decayME_,*outsideME_;  

};

//
// constants, enums and typedefs
//
using namespace std;
using namespace edm;

struct GlbSelectorStudy::MuonME {
  void bookHistograms(edm::Service<TFileService> fs, const std::string dirName)
  {
    TFileDirectory dir = fs->mkdir(dirName);
    TFileDirectory  StaDir = dir.mkdir( "Sta" );

    TH1F * h_pt = dir.make<TH1F>( "pt"  , "p_{t}", 100,  0., 100. );
    TH1F * h_pt2 = StaDir.make<TH1F>( "pt2"  , "p_{t}", 100,  0., 100. );
    h_pt3 = StaDir.make<TH1F>( "pt3"  , "p_{t}", 100,  0., 100. );

    //    TkDir = dir.mkdir( "Tk" );
    //    GlbDir = dir.mkdir( "Glb" );
  };
  void fill() {
    //
  };
  TH1F* h_pt3;
  //TFileDirectory StaDir;
  //  TH1 *dxy, *dz, *nHit, *chi2, *normChi2;
  //  TH1 *caloComp, *segComp;
};

//
// static data member definitions
//

//
// constructors and destructor
//
GlbSelectorStudy::GlbSelectorStudy(const edm::ParameterSet& iConfig):
  wantMotherBin(iConfig.getParameter<edm::ParameterSet>("IDconverttoBinNum"))
{
  theCategory = "GlbSelectorStudy";

   //now do what ever initialization is needed
  inTrackerBound_ = iConfig.getUntrackedParameter<bool>("inTrackerBound",false);
  inMuonBound_ = iConfig.getUntrackedParameter<bool>("inMuonBound",false);
  inCal_ = iConfig.getUntrackedParameter<bool>("inCal",false);
  simMuon_ =  iConfig.getUntrackedParameter<bool>("simMuon",false);

  //  glbMuLabel_ = iConfig.getParameter<edm::InputTag>("l3MuonLabel");
  //  trackingParticleLabel_ = iConfig.getParameter<edm::InputTag>("trackingParticleLabel");

  //  glbMuAssocLabel_ = iConfig.getParameter<std::string>("l3AssociatorName");

  // Labels for simulation and reconstruction tracks
  simLabel_  = iConfig.getParameter<InputTag>("simLabel" );
  trkMuLabel_ = iConfig.getParameter<InputTag>("trkMuLabel");
  staMuLabel_ = iConfig.getParameter<InputTag>("staMuLabel");
  glbMuLabel_ = iConfig.getParameter<InputTag>("glbMuLabel");
  muonLabel_ = iConfig.getParameter<InputTag>("muonLabel");

  // Labels for sim-reco association
  doAssoc_ = iConfig.getUntrackedParameter<bool>("doAssoc", true);
  trkMuAssocLabel_ = iConfig.getParameter<InputTag>("trkMuAssocLabel");
  staMuAssocLabel_ = iConfig.getParameter<InputTag>("staMuAssocLabel");
  glbMuAssocLabel_ = iConfig.getParameter<InputTag>("glbMuAssocLabel");


}


GlbSelectorStudy::~GlbSelectorStudy()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GlbSelectorStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;


  Handle<edm::SimTrackContainer> simTracks;
  iEvent.getByLabel("g4SimHits",simTracks);
  
  Handle<edm::SimVertexContainer> simVertexs;
  iEvent.getByLabel("g4SimHits",simVertexs);
  
  Handle<edm::HepMCProduct> hepmc;
  iEvent.getByType(hepmc);

  // Get TrackingParticles
  Handle<TrackingParticleCollection> simHandle;
  iEvent.getByLabel(simLabel_, simHandle);
  const TrackingParticleCollection simColl = *(simHandle.product());

  // Get Muon Tracks
  Handle<View<Track> > trkMuHandle;
  iEvent.getByLabel(trkMuLabel_, trkMuHandle);
  View<Track> trkMuColl = *(trkMuHandle.product());

  Handle<View<Track> > staMuHandle;
  iEvent.getByLabel(staMuLabel_, staMuHandle);
  View<Track> staMuColl = *(staMuHandle.product());

  Handle<View<Track> > glbMuHandle;
  iEvent.getByLabel(glbMuLabel_, glbMuHandle);
  View<Track> glbMuColl = *(glbMuHandle.product());

  // Get Muons
  Handle<View<Muon> > muonHandle;
  iEvent.getByLabel(muonLabel_, muonHandle);
  View<Muon> muonColl = *(muonHandle.product());

  // Get Association maps
  SimToRecoCollection simToTrkMuColl;
  SimToRecoCollection simToStaMuColl;
  SimToRecoCollection simToGlbMuColl;

  RecoToSimCollection trkMuToSimColl;
  RecoToSimCollection staMuToSimColl;
  RecoToSimCollection glbMuToSimColl;

  if ( doAssoc_ ) {
    // SimToReco associations
    simToTrkMuColl = trkMuAssociator_->associateSimToReco(trkMuHandle, simHandle, &iEvent);
    simToStaMuColl = staMuAssociator_->associateSimToReco(staMuHandle, simHandle, &iEvent);
    simToGlbMuColl = glbMuAssociator_->associateSimToReco(glbMuHandle, simHandle, &iEvent);
    
    // // RecoToSim associations
    trkMuToSimColl = trkMuAssociator_->associateRecoToSim(trkMuHandle, simHandle, &iEvent);
    staMuToSimColl = staMuAssociator_->associateRecoToSim(staMuHandle, simHandle, &iEvent);
    glbMuToSimColl = glbMuAssociator_->associateRecoToSim(glbMuHandle, simHandle, &iEvent);
  }
  else {
    // SimToReco associations
    Handle<SimToRecoCollection> simToTrkMuHandle;
    iEvent.getByLabel(trkMuAssocLabel_, simToTrkMuHandle);
    simToTrkMuColl = *(simToTrkMuHandle.product());
    
    Handle<SimToRecoCollection> simToStaMuHandle;
    iEvent.getByLabel(staMuAssocLabel_, simToStaMuHandle);
    simToStaMuColl = *(simToStaMuHandle.product());
    
    Handle<SimToRecoCollection> simToGlbMuHandle;
    iEvent.getByLabel(glbMuAssocLabel_, simToGlbMuHandle);
    simToGlbMuColl = *(simToGlbMuHandle.product());
    
    // RecoToSim associations
    Handle<RecoToSimCollection> trkMuToSimHandle;
    iEvent.getByLabel(trkMuAssocLabel_, trkMuToSimHandle);
    trkMuToSimColl = *(trkMuToSimHandle.product());
    
    Handle<RecoToSimCollection> staMuToSimHandle;
    iEvent.getByLabel(staMuAssocLabel_, staMuToSimHandle);
    staMuToSimColl = *(staMuToSimHandle.product());
    
    Handle<RecoToSimCollection> glbMuToSimHandle;
    iEvent.getByLabel(glbMuAssocLabel_, glbMuToSimHandle);
    glbMuToSimColl = *(glbMuToSimHandle.product());
  }
  
  int partype = 0;
  
  double decayR = 0;
  double decayZ = 0;
  
  bool inTracker = false;
  bool inMuon = false;
  
  bool simMuon = false;
  
  // Analyzer reco::Muon
  for(View<Muon>::const_iterator iMuon = muonColl.begin();
      iMuon != muonColl.end(); ++iMuon) {
 
    //const TrackRef glbTrack = (iMuon->isGlobalMuon()) ? iMuon->combinedMuon() : 0;
    
    if ( iMuon->isGlobalMuon() ) {
      const TrackRef glbTrack = iMuon->combinedMuon();
      //const RefToBase<Track> glbTrack = iMuon->combinedMuon();
      const RefToBase<Track> glbTrackRB(glbTrack);
    //}
    
      std::vector<std::pair<TrackingParticleRef,double> > tpRefV;
      if ( glbTrack.isAvailable() && trkMuToSimColl.find(glbTrackRB) != trkMuToSimColl.end() ) {
	tpRefV = trkMuToSimColl[glbTrackRB];

	const TrackingParticleRef & trp = tpRefV.begin()->first;

	int particle_ID = trp->pdgId();
	int myBin = wantMotherBin.GetBinNum(particle_ID);

	for(TrackingParticle::g4t_iterator isimtk = trp->g4Track_begin();isimtk!=trp->g4Track_end();isimtk++)
	  {//a
	
	    MotherSearch mother(&*isimtk, simTracks, simVertexs, hepmc);

		  if (mother.IsValid()){
		    if (mother.SimIsValid()){
		      //(*l3ParentID).push_back(mother.Sim_mother->type());
		      //(*l3MotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Sim_mother->type()));
		    }
		    else {
		      //(*l3ParentID).push_back(mother.Gen_mother->pdg_id());
		      //(*l3MotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id()));
		    }
		    //do it once per tracking particle once it succeed
		    break;
		  }
		  else{
		    edm::LogError(theCategory)<<"tricky muon from TrackingParticle.";
		  }
	  }//a
	

      }//find glbTrkRB
      
    }// isGlobal()

  }// loop over muon
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
GlbSelectorStudy::beginJob(const edm::EventSetup&)
{
edm::Service<TFileService> fs;

//TFileDirectory prompt = fs->mkdir( "prompt" );
//decay = fs->mkdir( "decay" );
//outside = fs->mkdir( "outside" );
 // TH1F * h_pt = outside.make<TH1F>( "pt"  , "p_{t}", 100,  0., 100. );

 promptME_ = new MuonME;
 // decayME_ = new MuonME;
 // outsideME_ = new MuonME;
 //
 promptME_->bookHistograms(fs,"prompt");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
GlbSelectorStudy::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(GlbSelectorStudy);
