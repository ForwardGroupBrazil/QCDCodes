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
// $Id: GlbSelectorStudy.cc,v 1.2 2009/04/02 16:26:32 aeverett Exp $
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

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include <TrackingTools/PatternTools/interface/TrajectoryMeasurement.h>
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
#include <TrackingTools/PatternTools/interface/Trajectory.h>
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include <DataFormats/TrackingRecHit/interface/TrackingRecHit.h>
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "TH1.h"
#include "TH2.h"
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
  virtual double kink(const reco::TrackRef& muon) const ;
  virtual void addTraj(const reco::Track& candIn);
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

  TrackTransformer *theTrackTransformer;

  //TFileDirectory * prompt; //, *promptSta, *promptTk, *promptGlb;
  //TFileDirectory decay;//, *decaySta, *decayTk, *decayGlb;
  //TFileDirectory outside;//, *outsideSta, *outsideTk, *outsideGlb;

  struct MuonME;
  MuonME *aME_,*bME_,*cME_;  
  MuonME *dME_,*eME_,*fME_, *gME_;  


};

//
// constants, enums and typedefs
//
using namespace std;
using namespace edm;
using namespace reco;

struct GlbSelectorStudy::MuonME {
  void bookHistograms(edm::Service<TFileService> fs, const std::string dirName)
  {
    TFileDirectory *dir = new TFileDirectory(fs->mkdir(dirName));
    TFileDirectory *StaDir = new TFileDirectory(dir->mkdir( "Sta" ));
    TFileDirectory *TkDir  = new TFileDirectory(dir->mkdir( "Trk" ));
    TFileDirectory *GlbDir = new TFileDirectory(dir->mkdir( "Glb" ));
    TFileDirectory *MuDir = new TFileDirectory(dir->mkdir( "Muon" ));
    /*
    vector<TFileDirectory*> dirV;
    dirV.push_back(StaDir);
    dirV.push_back(TkDir);
    dirV.push_back(GlbDir);

    vector<TFileDirectory*>::const_iterator dirIter;
    for(dirIter = dirV.begin(); dirIter!=dirV.end();++dirITer) {
      
    }
    */
    ht_dxy = TkDir->make<TH1F>("ht_dxy","d_{xy}",100,0.,1.);
    ht_dz = TkDir->make<TH1F>("ht_dz","d_{z}",100,0.,10.);
    ht_nHit = TkDir->make<TH1F>("ht_nHit","nHit",100,0.,100.);
    ht_chi2 = TkDir->make<TH1F>("ht_chi2","chi2",100,0.,100.);
    ht_nchi2 = TkDir->make<TH1F>("ht_nchi2","nchi2",100,0.,100.);

    hs_dxy = StaDir->make<TH1F>("hs_dxy","d_{xy}",100,0.,1.);
    hs_dz = StaDir->make<TH1F>("hs_dz","d_{z}",100,0.,10.);
    hs_nHit = StaDir->make<TH1F>("hs_nHit","nHit",100,0.,100.);
    hs_chi2 = StaDir->make<TH1F>("hs_chi2","chi2",100,0.,100.);
    hs_nchi2 = StaDir->make<TH1F>("hs_nchi2","nchi2",100,0.,100.);

    hg_dxy = GlbDir->make<TH1F>("hg_dxy","d_{xy}",100,0.,1.);
    hg_dz = GlbDir->make<TH1F>("hg_dz","d_{z}",100,0.,10.);
    hg_nHit = GlbDir->make<TH1F>("hg_nHit","nHit",100,0.,100.);
    hg_chi2 = GlbDir->make<TH1F>("hg_chi2","chi2",100,0.,100.);
    hg_nchi2 = GlbDir->make<TH1F>("hg_nchi2","nchi2",100,0.,100.);

    hm_hcal = MuDir->make<TH1F>("hm_hcal","E_{HCAL}",100,0.,10.);
    hm_ecal = MuDir->make<TH1F>("hm_ecal","E_{ECAL}",100,0.,10.);

    hs_NTrksEta_ = StaDir->make<TH1F>("hs_NTrksEta", "Number of reco tracks vs #eta", 50, -2.5, 2.5);
    hs_NTrksEta_St1_ = StaDir->make<TH1F>("hs_NTrksEta_St1", "Number of reco tracks vs #eta only in Station1", 50, -2.5, 2.5);
    hs_NTrksPt_ = StaDir->make<TH1F>("hs_NTrksPt", "Number of reco tracks vs p_{T}", 100, 0., 500.);
    hs_NTrksPt_St1_ = StaDir->make<TH1F>("hs_NTrksPt_St1", "Number of reco tracks vs p_{T} only in Station1", 100, 0., 500.);

    hg_NTrksEta_ = GlbDir->make<TH1F>("hg_NTrksEta", "Number of reco tracks vs #eta", 50, -2.5, 2.5);
    hg_NTrksEta_St1_ = GlbDir->make<TH1F>("hg_NTrksEta_St1", "Number of reco tracks vs #eta only in Station1", 50, -2.5, 2.5);
    hg_NTrksPt_ = GlbDir->make<TH1F>("hg_NTrksPt", "Number of reco tracks vs p_{T}", 100, 0., 500.);
    hg_NTrksPt_St1_ = GlbDir->make<TH1F>("hg_NTrksPt_St1", "Number of reco tracks vs p_{T} only in Station1", 100, 0., 500.);

    hm_nChamber = MuDir->make<TH1F>("hm_nChamber","nChamber",20,0,20);
    hm_nChamberMatch_no = MuDir->make<TH1F>("hm_nChamberMatch_no","nChamberMatch No Arbitration",20,0,20);
    hm_nChamberMatch_seg = MuDir->make<TH1F>("hm_nChamberMatch_seg","nChamberMatch Segment Arbitration",20,0,20);
    hm_nChamberMatch_segTrack = MuDir->make<TH1F>("hm_nChamberMatch_segTrack","nChamberMatch Segment and Track Arbitration",20,0,20);

    hm_caloComp = MuDir->make<TH1F>("hm_caloComp","caloComp",50,0.,1.);
    hm_segComp = MuDir->make<TH1F>("hm_segComp","segComp",50,0.,1.);
    hm_2DComp = MuDir->make<TH2F>("hm_2DComp","Seg/Calo Comp",50,0.,1.,50,0.,1.);

    hm_tm_sel = MuDir->make<TH1F>("hm_tm_sel","TM Selectors",11,-0.5,10.5);

    hm_outerPosition = MuDir->make<TH2F>("hm_outerPosition","OuterMost Position",1200,0.,1200.,1000,0.,1000.);
    hm_motherVtxPos = MuDir->make<TH2F>("hm_motherVtxPos","Mother Vtx Position",1201,-1.,1200.,1001,-1.,1000.);

    hg_kink = GlbDir->make<TH1F>("hg_kink","Kink",100,0.,10.);
    ht_kink = GlbDir->make<TH1F>("ht_kink","Kink",100,0.,10.);

    hm_tpType = MuDir->make<TH1F>("hm_tpType","TP Type",501,-0.5,500.5);
    hm_motherType = MuDir->make<TH1F>("hm_motherType","Mother Type",501,-0.5,500.5);

    hm_trackerMu = MuDir->make<TH1F>("hm_TM","isTrackerMuon",3,-1.5,1.5);

  };
  void fill(const reco::Muon& iMuon,const GlobalPoint pos) {
    //
    const TrackRef glbTrack = iMuon.combinedMuon();
    hg_dxy->Fill(glbTrack->dxy());
    hg_dz->Fill(glbTrack->dz());
    hg_nHit->Fill(glbTrack->numberOfValidHits());
    hg_chi2->Fill(glbTrack->chi2());
    hg_nchi2->Fill(glbTrack->normalizedChi2());

    hm_trackerMu->Fill(iMuon.isTrackerMuon());

    hg_NTrksEta_->Fill(glbTrack->eta());
    hg_NTrksPt_->Fill(glbTrack->pt());
    int station = 0;
    DetId id(glbTrack->outerDetId());
    if ( id.subdetId() == MuonSubdetId::DT ) {
      DTChamberId did(id.rawId());
      station = did.station();
    }  else if ( id.subdetId() == MuonSubdetId::CSC ) {
      CSCDetId did(id.rawId());
      station = did.station();
    }   else if ( id.subdetId() == MuonSubdetId::RPC ) {
      RPCDetId rpcid(id.rawId());
      station = rpcid.station();
    }
    
    if(station == 1) hg_NTrksEta_St1_->Fill(glbTrack->eta());
    if(station == 1) hg_NTrksPt_St1_->Fill(glbTrack->pt());

    const TrackRef staTrack = iMuon.standAloneMuon();
    hs_dxy->Fill(staTrack->dxy());
    hs_dz->Fill(staTrack->dz());
    hs_nHit->Fill(staTrack->numberOfValidHits());
    hs_chi2->Fill(staTrack->chi2());
    hs_nchi2->Fill(staTrack->normalizedChi2());

    hs_NTrksEta_->Fill(staTrack->eta());
    hs_NTrksPt_->Fill(staTrack->pt());
    station = 0;
    DetId id2(staTrack->outerDetId());
    if ( id2.subdetId() == MuonSubdetId::DT ) {
      DTChamberId did(id2.rawId());
      station = did.station();
    }  else if ( id2.subdetId() == MuonSubdetId::CSC ) {
      CSCDetId did(id2.rawId());
      station = did.station();
    }   else if ( id2.subdetId() == MuonSubdetId::RPC ) {
      RPCDetId rpcid(id2.rawId());
      station = rpcid.station();
    }
    
    if(station == 1) hs_NTrksEta_St1_->Fill(staTrack->eta());
    if(station == 1) hs_NTrksPt_St1_->Fill(staTrack->pt());
    
    const TrackRef trkTrack = iMuon.track();
    ht_dxy->Fill(trkTrack->dxy());
    ht_dz->Fill(trkTrack->dz());
    ht_nHit->Fill(trkTrack->numberOfValidHits());
    ht_chi2->Fill(trkTrack->chi2());
    ht_nchi2->Fill(trkTrack->normalizedChi2());



    if(iMuon.isEnergyValid()) hm_hcal->Fill(iMuon.calEnergy().had);
    if(iMuon.isEnergyValid()) hm_ecal->Fill(iMuon.calEnergy().em);

    hm_nChamber->Fill(iMuon.numberOfChambers());
    hm_nChamberMatch_no->Fill(iMuon.numberOfMatches(reco::Muon::NoArbitration));
    hm_nChamberMatch_seg->Fill(iMuon.numberOfMatches(reco::Muon::SegmentArbitration));
    hm_nChamberMatch_segTrack->Fill(iMuon.numberOfMatches(reco::Muon::SegmentAndTrackArbitration));

    if(iMuon.isCaloCompatibilityValid()) hm_caloComp->Fill(iMuon.caloCompatibility());
    hm_segComp->Fill(iMuon.segmentCompatibility());
    if(iMuon.isCaloCompatibilityValid()) 
      hm_2DComp->Fill(iMuon.caloCompatibility(),iMuon.segmentCompatibility());
    else
      hm_2DComp->Fill(0.,iMuon.segmentCompatibility());

    
    hm_tm_sel->Fill(0);
    if(iMuon.isGood(reco::Muon::TMLastStationLoose)) hm_tm_sel->Fill(1);
    if(iMuon.isGood(reco::Muon::TMLastStationTight)) hm_tm_sel->Fill(2);
    if(iMuon.isGood(reco::Muon::TM2DCompatibilityLoose)) hm_tm_sel->Fill(3);
    if(iMuon.isGood(reco::Muon::TM2DCompatibilityTight)) hm_tm_sel->Fill(4);
    if(iMuon.isGood(reco::Muon::TMOneStationLoose)) hm_tm_sel->Fill(5);
    if(iMuon.isGood(reco::Muon::TMOneStationTight)) hm_tm_sel->Fill(6);
    if(iMuon.isGood(reco::Muon::AllTrackerMuons)) hm_tm_sel->Fill(7);
    if(iMuon.isGood(reco::Muon::TrackerMuonArbitrated)) hm_tm_sel->Fill(8);
    if(iMuon.isGood(reco::Muon::TMLastStationOptimizedLowPtLoose)) hm_tm_sel->Fill(7);
    if(iMuon.isGood(reco::Muon::TMLastStationOptimizedLowPtTight)) hm_tm_sel->Fill(8);

    float outerX = glbTrack->outerX();
    float outerY = glbTrack->outerY();
    float outerZ = glbTrack->outerZ();
    float outerR = sqrt( outerX*outerX +outerY*outerY);
    hm_outerPosition->Fill(outerZ,outerR);
    hm_motherVtxPos->Fill(abs(pos.z()),pos.perp());


  };


  TH1F *ht_dxy, *ht_dz, *ht_nHit, *ht_chi2, *ht_nchi2;
  TH1F *hs_dxy, *hs_dz, *hs_nHit, *hs_chi2, *hs_nchi2;
  TH1F *hg_dxy, *hg_dz, *hg_nHit, *hg_chi2, *hg_nchi2;

  TH1F *hm_hcal, *hm_ecal;

  TH1F *hs_NTrksEta_, *hs_NTrksEta_St1_,  *hs_NTrksPt_,  *hs_NTrksPt_St1_;
  TH1F *hg_NTrksEta_, *hg_NTrksEta_St1_,  *hg_NTrksPt_,  *hg_NTrksPt_St1_;

  TH1F *hm_nChamber;
  TH1F *hm_nChamberMatch_no, *hm_nChamberMatch_seg, *hm_nChamberMatch_segTrack;
  TH1F *hm_caloComp, *hm_segComp;

  TH1F *hm_tm_sel;

  TH2F *hm_outerPosition;
  TH2F *hm_motherVtxPos;
  TH2F *hm_2DComp;

  TH1F *hg_kink;
  TH1F *ht_kink;

  TH1F *hm_tpType, *hm_motherType, *hm_trackerMu;
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

  ParameterSet trackTransformerPSet = iConfig.getParameter<ParameterSet>("TrackTransformer");
  //trackTransformerPSet.addParameter<string>("Propagator",TransformerOutPropagator);
  //theTrackTransformer = new TrackTransformer(trackTransformerPSet);

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
  LogDebug(theCategory);

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
  


  LogDebug(theCategory)<<"MuonColl size " << muonColl.size();  
  // Analyzer reco::Muon
  for(View<Muon>::const_iterator iMuon = muonColl.begin();
      iMuon != muonColl.end(); ++iMuon) {

    int muClass = 0;
    int tpType = 0;
    int motherType = 0;
    bool simMuon = false;
    
    double decayR = 0;
    double decayZ = 0;
    
    bool inTracker = false;
    bool inMuon = false;
    bool validMother = false;
    bool validSimMother = false;
    bool validGenMother = false;
        
    double x = -1.0;
    double y = -1.0;
    double z = -1.0; 

    LogTrace(theCategory)<<"Looking at a new muon.....";
    if ( iMuon->isGlobalMuon() ) {      
      LogTrace(theCategory)<<"which isGlobalMuon";

      const TrackRef glbTrack = iMuon->combinedMuon();
      const RefToBase<Track> glbTrackRB(glbTrack);
      
      LogTrace(theCategory)<<"that is available " <<  glbTrack.isAvailable();
      
      std::vector<std::pair<TrackingParticleRef,double> > tpRefV;
      if ( glbTrack.isAvailable() && glbMuToSimColl.find(glbTrackRB) != glbMuToSimColl.end() ) {//get TP
	tpRefV = glbMuToSimColl[glbTrackRB];

	LogTrace(theCategory)<<"Found tpRefV of size " << tpRefV.size();

	const TrackingParticleRef & trp = tpRefV.begin()->first;

	int particle_ID = trp->pdgId();
	tpType = particle_ID;
	int myBin = wantMotherBin.GetBinNum(particle_ID);

	LogTrace(theCategory)<<"with the leading TP ID " << particle_ID;
	
	for(TrackingParticle::g4t_iterator isimtk = trp->g4Track_begin();isimtk!=trp->g4Track_end();isimtk++) {//loop over sim in TP
	  LogTrace(theCategory)<<"... now going to look for mother ....";
	  MotherSearch mother(&*isimtk, simTracks, simVertexs, hepmc);

	  validMother = mother.IsValid();	  
	  if (validMother){
	    validSimMother = mother.SimIsValid();
	    validGenMother = mother.GenIsValid();
	    if (validSimMother){
	      LogTrace(theCategory)<<"motherFromSim " << mother.Sim_mother->type();
	      LogTrace(theCategory)<<"     vertex " << mother.Sim_vertex->position().Rho() << " " <<  mother.Sim_vertex->position().z();
	      x = mother.Sim_vertex->position().x();
	      y = mother.Sim_vertex->position().y();
	      z = mother.Sim_vertex->position().z();
	      motherType = mother.Sim_mother->type();
	      //(*l3ParentID).push_back(mother.Sim_mother->type());
	      //(*l3MotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Sim_mother->type()));
	    } else {
	      LogTrace(theCategory)<<"motherFromGen " << mother.Gen_mother->pdg_id();
	      LogTrace(theCategory)<<"     vertex " << mother.Gen_vertex->position().perp() << " " <<  mother.Gen_vertex->position().z();
	      x = 0.1 * mother.Gen_vertex->position().x();
	      y = 0.1 * mother.Gen_vertex->position().y();
	      z = 0.1 * mother.Gen_vertex->position().z();
	      motherType = mother.Gen_mother->pdg_id();
	      //(*l3ParentID).push_back(mother.Gen_mother->pdg_id());
	      //(*l3MotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id()));
	    }
	    //do it once per tracking particle once it succeed
	    break;
	  } else{
	    //this is a "prompt" TrackingParticle
	    edm::LogError(theCategory)<<"tricky muon from TrackingParticle.";
	  }
	}//loop over sim in TP
      }//get TP
      GlobalPoint pos(x,y,z);

      if(TrackerBounds::isInside(pos)) inTracker = true;
      if(MuonBounds::isInside(pos)) inMuon = true;
      LogTrace(theCategory) << "***TP is inMuon " << inMuon << " and inTracker " << inTracker;

      if(!validMother) {
	if(abs(tpType)==13) muClass = 1; //prompt mu
	else muClass = 2; //prompt other (punch-through)
      } else {
	if(abs(tpType)==13 && abs(motherType)!=13) {
	  if(inTracker) muClass = 3; //decay to muon
	  else muClass = 4; //shower to muon
	} else if (abs(tpType)==13 && abs(motherType)==13) {
	  muClass = 5; //muon decay to muon
	} else if (abs(tpType)!=13) {
	  if(inTracker) muClass = 6; //othr (prompt other? : punch-through?)
	  else muClass = 7;
	}
      }

      LogTrace(theCategory)<<"MuClass " << muClass;

      /*
      if(muClass==1) aME_->fill(*iMuon,pos);
      if(muClass==2) bME_->fill(*iMuon,pos);
      if(muClass==3) cME_->fill(*iMuon,pos);
      if(muClass==4) dME_->fill(*iMuon,pos);
      if(muClass==5) eME_->fill(*iMuon,pos);
      if(muClass==6) fME_->fill(*iMuon,pos);
      */

      MuonME * thisME = 0;

      if(muClass==1) thisME = aME_;
      if(muClass==2) thisME = bME_;
      if(muClass==3) thisME = cME_;
      if(muClass==4) thisME = dME_;
      if(muClass==5) thisME = eME_;
      if(muClass==6) thisME = fME_;
      if(muClass==7) thisME = gME_;

      thisME->fill(*iMuon,pos);
      //thisME->hg_kink->Fill(kink(iMuon->combinedMuon()));
      //thisME->ht_kink->Fill(kink(iMuon->innerTrack()));
      thisME->hm_tpType->Fill(abs(tpType));
      thisME->hm_motherType->Fill(abs(motherType));

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

 aME_ = new MuonME;
 bME_ = new MuonME;
 cME_ = new MuonME;
 dME_ = new MuonME;
 eME_ = new MuonME;
 fME_ = new MuonME;
 gME_ = new MuonME;
 // decayME_ = new MuonME;
 // outsideME_ = new MuonME;
 //
 aME_->bookHistograms(fs,"A");
 bME_->bookHistograms(fs,"B");
 cME_->bookHistograms(fs,"C");
 dME_->bookHistograms(fs,"D");
 eME_->bookHistograms(fs,"E");
 fME_->bookHistograms(fs,"F");
 gME_->bookHistograms(fs,"G");



}

// ------------ method called once each job just after ending the event loop  ------------
void 
GlbSelectorStudy::endJob() {
}


//
// kink finder
//
void GlbSelectorStudy::addTraj(const reco::Track& candIn) {

  //if ( candIn.first == 0 ) {
  
  //Trajectory *returnTrajectory = 0;
  typedef std::vector<Trajectory> TC;
  
  TC staTrajs = theTrackTransformer->transform(candIn);

  if (staTrajs.empty()) {
    LogTrace(theCategory) << "Transformer: Add Traj failed!";
    //candIn = TrackCand(0,candIn.second); 
  } else {
    //return staTrajs.front();
    //Trajectory * tmpTrajectory = new Trajectory(staTrajs.front());
    //tmpTrajectory->setSeedRef(candIn.second->seedRef());
    //candIn = TrackCand(tmpTrajectory,candIn.second);
  }
  //}
}

double GlbSelectorStudy::kink(const reco::TrackRef& muonIn) const {
    
    typedef std::vector<Trajectory> TC;
    
    TC staTrajs = theTrackTransformer->transform(muonIn);
    
    if (staTrajs.empty()) {
      LogTrace(theCategory) << "Transformer: Add Traj failed!";
      return 0.;
    } 
      Trajectory muon(staTrajs.front());
    
    
    double result = 0.0;
    
    typedef TransientTrackingRecHit::ConstRecHitPointer 	ConstRecHitPointer;
    typedef ConstRecHitPointer RecHit;
    typedef vector<TrajectoryMeasurement>::const_iterator TMI;
    
    vector<TrajectoryMeasurement> meas = muon.measurements();
    
    for ( TMI m = meas.begin(); m != meas.end(); m++ ) {
      RecHit rhit = (*m).recHit();
      bool ok = false;
      if ( rhit->isValid() ) {
	//const GeomDetType& type = rhit->det()->detUnits().front()->type();
	//if ( type.module() == pixel || type.module() == silicon ) ok = true;
	
	if(DetId::Tracker == rhit->geographicalId().det()) ok = true;
	
      }
      if ( !ok ) continue;
      
      const TrajectoryStateOnSurface& tsos = (*m).predictedState();
      if ( tsos.isValid() ) {
	double phi1 = tsos.globalPosition().phi();
	if ( phi1 < 0 ) phi1 = 2*M_PI + phi1;
	
	double phi2 = rhit->globalPosition().phi();
	if ( phi2 < 0 ) phi2 = 2*M_PI + phi2;
	
	double diff = fabs(phi1 - phi2);
	if ( diff > M_PI ) diff = 2*M_PI - diff;
	
	//GlobalPoint hitPos = rhit->det()->toGlobal(rhit->localPosition());
	//GlobalError hitErr = rhit->det()->toGlobal(rhit->localPositionError());
	GlobalPoint hitPos = rhit->globalPosition();
	GlobalError hitErr = rhit->globalPositionError();
	
	double error = hitErr.phierr(hitPos);  // error squared
	
	double s = ( error > 0.0 ) ? (diff*diff)/error : (diff*diff);
	result += s;
	
      }
    }
    
    return result;
    
  }



//define this as a plug-in
DEFINE_FWK_MODULE(GlbSelectorStudy);
