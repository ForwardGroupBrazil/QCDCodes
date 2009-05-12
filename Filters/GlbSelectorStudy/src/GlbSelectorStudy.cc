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
// $Id: GlbSelectorStudy.cc,v 1.4 2009/04/15 19:24:55 aeverett Exp $
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

#include "PhysicsTools/RecoAlgos/interface/TrackingParticleSelector.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackHistory/interface/TrackClassifier.h"
#include "SimTracker/TrackHistory/interface/TrackCategories.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

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
  
  static const unsigned int primaryMuon       =  1<<1;
  static const unsigned int siliconMuon       =  1<<2;
  static const unsigned int calConversionMuon =  1<<3;
  static const unsigned int otherMuon         =  1<<4;
  static const unsigned int punchThrough      =  1<<5;
  
  //void setType( unsigned int type ) { type_ = type; }
  //unsigned int type() const { return type_; }
  
  bool isPrimaryMuon(unsigned int inType_)     const { return  inType_ & primaryMuon; }
  bool isSiliconMuon(unsigned int inType_)    const { return  inType_ & siliconMuon; }
  bool isCalConversionMuon(unsigned int inType_) const { return  inType_ & calConversionMuon; }
  bool isOtherMuon(unsigned int inType_) const { return  inType_ & otherMuon; }
  bool isPunchThrough(unsigned int inType_) const { return  inType_ & punchThrough; }
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual double kink(const reco::TrackRef& muon) const ;
  virtual void addTraj(const reco::Track& candIn);
  virtual void printTruth(TrackingParticleRef & simRef, int i);
  virtual unsigned int getBit(const TrackingParticleRef&) const;
  virtual std::string particleString(int pdgId) const;
  virtual std::string vertexString(TrackingParticleRefVector in, TrackingParticleRefVector out ) const;
  // ----------member data ---------------------------
  //unsigned int type_;
  
  bool inTrackerBound_;
  bool inMuonBound_;
  bool inCal_;
  bool simMuon_;
  
  bool doAssoc_;

    edm::ESHandle<ParticleDataTable> pdt_;
  
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

  TrackingParticleSelector tpSelector_primary;
  TrackingParticleSelector tpSelector_silicon;
  TrackingParticleSelector tpSelector_calConversion;
  TrackClassifier classifier_;

};

//
// constants, enums and typedefs
//
typedef std::vector<TrackingVertex>                TrackingVertexCollection;
typedef edm::Ref<TrackingVertexCollection>         TrackingVertexRef;
typedef edm::RefVector<TrackingVertexCollection>   TrackingVertexRefVector;
typedef TrackingVertexRefVector::iterator          tv_iterator;

typedef edm::AssociationMap<edm::OneToValue<TrackingParticleCollection, unsigned int, unsigned int > > tpBitMap;
typedef edm::AssociationMap<edm::OneToValue<reco::TrackCollection, unsigned int, unsigned int > > recoBitMap;
typedef edm::AssociationMap<edm::OneToOne<reco::TrackCollection, TrackingParticleCollection, unsigned int > > recoTpMap;

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
    ht_nchi2Prob = TkDir->make<TH1F>("ht_nchi2Prob","nchi2 Probability",100,0.,100.);

    hs_dxy = StaDir->make<TH1F>("hs_dxy","d_{xy}",100,0.,1.);
    hs_dz = StaDir->make<TH1F>("hs_dz","d_{z}",100,0.,10.);
    hs_nHit = StaDir->make<TH1F>("hs_nHit","nHit",100,0.,100.);
    hs_chi2 = StaDir->make<TH1F>("hs_chi2","chi2",100,0.,100.);
    hs_nchi2 = StaDir->make<TH1F>("hs_nchi2","nchi2",100,0.,100.);
    hs_nchi2Prob = StaDir->make<TH1F>("hs_nchi2Prob","nchi2 Probability",100,0.,100.);

    hg_dxy = GlbDir->make<TH1F>("hg_dxy","d_{xy}",100,0.,1.);
    hg_dz = GlbDir->make<TH1F>("hg_dz","d_{z}",100,0.,10.);
    hg_nHit = GlbDir->make<TH1F>("hg_nHit","nHit",100,0.,100.);
    hg_chi2 = GlbDir->make<TH1F>("hg_chi2","chi2",100,0.,100.);
    hg_nchi2 = GlbDir->make<TH1F>("hg_nchi2","nchi2",100,0.,100.);
    hg_nchi2Prob = GlbDir->make<TH1F>("hg_nchi2Prob","nchi2 Probability",100,0.,100.);

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

    hm_NSt_tot_ =MuDir->make<TH1F>("hm_NSt_tot","nStation total",11,-0.5,10.5);
    hm_NSt_rpc_ =MuDir->make<TH1F>("hm_NSt_rpc","nStation RPC",11,-0.5,10.5);
    hm_NSt_csc_ =MuDir->make<TH1F>("hm_NSt_csc","nStation CSC",11,-0.5,10.5);
    hm_NSt_dt_ =MuDir->make<TH1F>("hm_NSt_dt","nStation DT",11,-0.5,10.5);

    hm_noDepth1_ = MuDir->make<TH1F>("hm_noDepth1","no depth (1)",3,-1.5,1.5);
    hm_noDepth2_ = MuDir->make<TH1F>("hm_noDepth2","no depth (2)",3,-1.5,1.5);

    hm_outerPosition = MuDir->make<TH2F>("hm_outerPosition","OuterMost Position",1200,0.,1200.,1000,0.,1000.);
    hm_motherVtxPos = MuDir->make<TH2F>("hm_motherVtxPos","Mother Vtx Position",1201,-1.,1200.,1001,-1.,1000.);
    hm_motherDecayPos = MuDir->make<TH2F>("hm_motherDecayPos","Mother Decay Position",1201,-1.,1200.,1001,-1.,1000.);
    hm_motherMuVtxPos = MuDir->make<TH2F>("hm_motherMuVtxPos","MotherMu Vtx Position",1201,-1.,1200.,1001,-1.,1000.);
    hm_motherMuDecayPos = MuDir->make<TH2F>("hm_motherMuDecayPos","MotherMu Decay Position",1201,-1.,1200.,1001,-1.,1000.);

    hg_kink = GlbDir->make<TH1F>("hg_kink","Kink",100,0.,10.);
    ht_kink = TkDir->make<TH1F>("ht_kink","Kink",100,0.,10.);

    hm_tpType = MuDir->make<TH1F>("hm_tpType","TP Type",501,-0.5,500.5);
    hm_motherType = MuDir->make<TH1F>("hm_motherType","Mother Type",501,-0.5,500.5);

    hm_trackerMu = MuDir->make<TH1F>("hm_TM","isTrackerMuon",3,-1.5,1.5);

  };
  //void fill(const reco::Muon& iMuon,const GlobalPoint pos,const GlobalPoint decayPos) {
  void fill(const reco::Muon& iMuon,const TrackingParticleRef& pos, const TrackingParticleRef& muTp) {
    //
    hm_tpType->Fill(pos->pdgId());

    const TrackRef glbTrack = iMuon.combinedMuon();
    hg_dxy->Fill(glbTrack->dxy());
    hg_dz->Fill(glbTrack->dz());
    hg_nHit->Fill(glbTrack->numberOfValidHits());
    hg_chi2->Fill(glbTrack->chi2());
    hg_nchi2->Fill(glbTrack->normalizedChi2());
    hg_nchi2Prob->Fill(-LnChiSquaredProbability(glbTrack->chi2(), glbTrack->ndof()));


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
    hs_nchi2Prob->Fill(-LnChiSquaredProbability(staTrack->chi2(), staTrack->ndof()));

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

    std::pair<std::pair<int,int>,std::pair<int,int> > nStat = countStations(iMuon.combinedMuon());

    hm_NSt_tot_->Fill(nStat.first.first);
    hm_NSt_rpc_->Fill(nStat.first.second);
    hm_NSt_dt_->Fill(nStat.second.first);
    hm_NSt_csc_->Fill(nStat.second.second);

    bool noDepth = (nStat.first.second == 1 && (nStat.second.first == 1 || nStat.second.second == 1) ) ? 1 : 0;

    int nCSCSeg[4];
    int nDTSeg[4];
    int nRPCSeg[4];
    for(int station = 0; station < 4; ++station) {
      nCSCSeg[station] = iMuon.numberOfSegments(station+1, MuonSubdetId::CSC, Muon::NoArbitration);
      nDTSeg[station] = iMuon.numberOfSegments(station+1, MuonSubdetId::DT, Muon::NoArbitration);
      nRPCSeg[station] = iMuon.numberOfSegments(station+1, MuonSubdetId::RPC, Muon::NoArbitration);
    }

    bool noDepth2 =  ( (nCSCSeg[0] > 1 && nCSCSeg[1] == 0 && nCSCSeg[2] == 0 && nCSCSeg[3] == 0 && nRPCSeg[0] >= 1) ||
	 (nDTSeg[0] > 1 && nDTSeg[1] == 0 && nDTSeg[2] == 0 && nDTSeg[3] == 0 && nRPCSeg[0] >= 1) ) ? 1 : 0;

    hm_noDepth1_->Fill(noDepth);
    hm_noDepth2_->Fill(noDepth2);


    
    const TrackRef trkTrack = iMuon.track();
    ht_dxy->Fill(trkTrack->dxy());
    ht_dz->Fill(trkTrack->dz());
    ht_nHit->Fill(trkTrack->numberOfValidHits());
    ht_chi2->Fill(trkTrack->chi2());
    ht_nchi2->Fill(trkTrack->normalizedChi2());
    ht_nchi2Prob->Fill(-LnChiSquaredProbability(trkTrack->chi2(), trkTrack->ndof()));


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

    if(pos.isAvailable()) {
      tv_iterator dv = pos->decayVertices_begin();
      float decayx = (*dv)->position().x();
      float decayy = (*dv)->position().y();
      float decayz = (*dv)->position().z();
      
      hm_motherVtxPos->Fill(abs(pos->vertex().z()),sqrt(pos->vertex().perp2()));
      hm_motherDecayPos->Fill(abs(decayz),sqrt(decayx*decayx + decayy*decayy) );

    }

    if(muTp.isAvailable()) {
      tv_iterator dv = muTp->decayVertices_begin();
      float decayx = (*dv)->position().x();
      float decayy = (*dv)->position().y();
      float decayz = (*dv)->position().z();
      
      hm_motherMuVtxPos->Fill(abs(muTp->vertex().z()),sqrt(muTp->vertex().perp2()));
      hm_motherMuDecayPos->Fill(abs(decayz),sqrt(decayx*decayx + decayy*decayy) );

    }

  };

std::pair< std::pair<int,int>,std::pair<int,int> > countStations(const reco::TrackRef& track)
{
  bool dt_1 = false;
  bool dt_2 = false;
  bool dt_3 = false;
  bool dt_4 = false;
  bool csc_1 = false;
  bool csc_2 = false;
  bool csc_3 = false;
  bool csc_4 = false;
  bool rpc_1 = false;
  bool rpc_2 = false;
  bool rpc_3 = false;
  bool rpc_4 = false;

  for (trackingRecHit_iterator iall = track->recHitsBegin(); iall != track->recHitsEnd(); ++iall) {
    
    //    if( (*iall)->det()->geographicalId().det()==2 && (*iall)->det()->geographicalId().subdetId()==3) continue;

    if((*iall)->isValid()){
      int DT_station = 0;
      int CSC_station = 0;
      int RPC_station = 0;
      
      switch((*iall)->geographicalId().det()) {      
      case DetId::Tracker:
        break;
      case DetId::Muon:      
        switch((*iall)->geographicalId().subdetId()) {  
        case MuonSubdetId::DT:
          {
          DTChamberId detId_dt((*iall)->geographicalId().rawId());
          DT_station = detId_dt.station();
          break;
          }
        case MuonSubdetId::CSC:
          {
	    CSCDetId detId_csc((*iall)->geographicalId().rawId());
	    CSC_station = detId_csc.station();
	    break;
          }
	case MuonSubdetId::RPC:
	  {
	    RPCDetId rpcid((*iall)->geographicalId().rawId());
	    RPC_station = rpcid.station();	    
	    break;
	  }
        default: ;
        }
        break;
      default: ;
      }
      
      if(DT_station == 1) dt_1 = true;
      if(DT_station == 2) dt_2 = true;
      if(DT_station == 3) dt_3 = true;
      if(DT_station == 4) dt_4 = true;
      
      if(CSC_station == 1) csc_1 = true;
      if(CSC_station == 2) csc_2 = true;
      if(CSC_station == 3) csc_3 = true;
      if(CSC_station == 4) csc_4 = true;

      if(RPC_station == 1) rpc_1 = true;
      if(RPC_station == 2) rpc_2 = true;
      if(RPC_station == 3) rpc_3 = true;
      if(RPC_station == 4) rpc_4 = true;
      
    }
  }
  int last = 0;
  int n_dt = 0;
  if (dt_1) {n_dt++;last=1;}
  if (dt_2) {n_dt++;last=2;}
  if (dt_3) {n_dt++;last=3;}
  if (dt_4) {n_dt++;last=4;}

  int n_csc = 0;
  if (csc_1) {n_csc++;last=5;}
  if (csc_2) {n_csc++;last=6;}
  if (csc_3) {n_csc++;last=7;}
  if (csc_4) {n_csc++;last=8;}

  int n_rpc = 0;
  if (rpc_1) {n_rpc++;last=9;}
  if (rpc_2) {n_rpc++;last=10;}
  if (rpc_3) {n_rpc++;last=11;}
  if (rpc_4) {n_rpc++;last=12;}
  
  int nTot = (n_dt+n_csc+n_rpc);

  std::pair<int,int> p1(nTot,n_rpc);
  std::pair<int,int> p2(n_dt,n_csc);
 
  return std::pair<std::pair<int,int>,std::pair<int,int> >(p1,p2);
};


  TH1F *ht_dxy, *ht_dz, *ht_nHit, *ht_chi2, *ht_nchi2, *ht_nchi2Prob;
  TH1F *hs_dxy, *hs_dz, *hs_nHit, *hs_chi2, *hs_nchi2, *hs_nchi2Prob;
  TH1F *hg_dxy, *hg_dz, *hg_nHit, *hg_chi2, *hg_nchi2, *hg_nchi2Prob;

  TH1F *hm_hcal, *hm_ecal;

  TH1F *hs_NTrksEta_, *hs_NTrksEta_St1_,  *hs_NTrksPt_,  *hs_NTrksPt_St1_;
  TH1F *hg_NTrksEta_, *hg_NTrksEta_St1_,  *hg_NTrksPt_,  *hg_NTrksPt_St1_;

  TH1F *hm_nChamber;
  TH1F *hm_nChamberMatch_no, *hm_nChamberMatch_seg, *hm_nChamberMatch_segTrack;
  TH1F *hm_caloComp, *hm_segComp;

  TH1F *hm_tm_sel;

  TH2F *hm_outerPosition;
  TH2F *hm_motherVtxPos,*hm_motherDecayPos;
  TH2F *hm_motherMuVtxPos,*hm_motherMuDecayPos;
  TH2F *hm_2DComp;

  TH1F *hg_kink;
  TH1F *ht_kink;

  TH1F *hm_tpType, *hm_motherType, *hm_trackerMu;

  TH1F *hm_NSt_tot_, *hm_NSt_rpc_, *hm_NSt_csc_, *hm_NSt_dt_;
  TH1F *hm_noDepth1_, *hm_noDepth2_;
};

//
// static data member definitions
//

//
// constructors and destructor
//
GlbSelectorStudy::GlbSelectorStudy(const edm::ParameterSet& iConfig):
  wantMotherBin(iConfig.getParameter<edm::ParameterSet>("IDconverttoBinNum")),
  classifier_(iConfig)
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
  theTrackTransformer = new TrackTransformer(trackTransformerPSet);

  ParameterSet tpset_primary = iConfig.getParameter<ParameterSet>("tpSelector_primary");
  tpSelector_primary = TrackingParticleSelector(					       tpset_primary.getParameter<double>("ptMin"),
											       tpset_primary.getParameter<double>("minRapidity"),
											       tpset_primary.getParameter<double>("maxRapidity"),
											       tpset_primary.getParameter<double>("tip"),
											       tpset_primary.getParameter<double>("lip"),
											       tpset_primary.getParameter<int>("minHit"),
											       tpset_primary.getParameter<bool>("signalOnly"),
											       tpset_primary.getParameter<bool>("chargedOnly"),
											       tpset_primary.getParameter<std::vector<int> >("pdgId"));
  
  ParameterSet tpset_silicon = iConfig.getParameter<ParameterSet>("tpSelector_silicon");
  tpSelector_silicon = TrackingParticleSelector(tpset_silicon.getParameter<double>("ptMin"),
                                         tpset_silicon.getParameter<double>("minRapidity"),
                                         tpset_silicon.getParameter<double>("maxRapidity"),
                                         tpset_silicon.getParameter<double>("tip"),
                                         tpset_silicon.getParameter<double>("lip"),
                                         tpset_silicon.getParameter<int>("minHit"), 
                                         tpset_silicon.getParameter<bool>("signalOnly"),
                                         tpset_silicon.getParameter<bool>("chargedOnly"),
                                         tpset_silicon.getParameter<std::vector<int> >("pdgId"));

  ParameterSet tpset_calConversion = iConfig.getParameter<ParameterSet>("tpSelector_calConversion");
  tpSelector_calConversion = TrackingParticleSelector(tpset_calConversion.getParameter<double>("ptMin"),
                                         tpset_calConversion.getParameter<double>("minRapidity"),
                                         tpset_calConversion.getParameter<double>("maxRapidity"),
                                         tpset_calConversion.getParameter<double>("tip"),
                                         tpset_calConversion.getParameter<double>("lip"),
                                         tpset_calConversion.getParameter<int>("minHit"),
                                         tpset_calConversion.getParameter<bool>("signalOnly"),
                                         tpset_calConversion.getParameter<bool>("chargedOnly"),
                                         tpset_calConversion.getParameter<std::vector<int> >("pdgId"));
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

  classifier_.newEvent(iEvent, iSetup);



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

  ////////////////////////////////



  ////typedef  std::map< reco::TrackCollection, unsigned int > myMap1;
  //typedef  std::map< unsigned int, unsigned int > myMap2;
  
  //typedef edm::AssociationMap<
  //edm::OneToOne<reco::TrackCollection,
  //reco::TrackCollection> > myOtherMap;
  
  
  tpBitMap   tpBitMap_;
  recoBitMap recoBitMap_;
  recoTpMap recoTpMap_;
  //myMap2 testMap2;
  //myOtherMap testMap3;
  
    //typedef std::vector<int> CKey;
    //typedef std::vector<double> CVal;
    //typedef edm::AssociationMap<edm::OneToOne<CKey, CVal, unsigned char> > Assoc;


  ///////////////////////
  LogDebug(theCategory)<<"SimColl size " << simColl.size();
  const TrackingParticleCollection::size_type nSim = simColl.size();
  for(TrackingParticleCollection::size_type i=0; i<nSim; i++) {
    TrackingParticleRef simRef(simHandle, i);
    const TrackingParticle* simTP = simRef.get();
    
    unsigned int thisBit = getBit(simRef);

    if ( ! (isPrimaryMuon(thisBit) || isSiliconMuon(thisBit) || isCalConversionMuon(thisBit) ) ) continue;

    printTruth(simRef, i);
    
    std::vector<std::pair<RefToBase<Track>, double> > rt;
    
    if(simToGlbMuColl.find(simRef) != simToGlbMuColl.end()){
      //rt = (std::vector<std::pair<RefToBase<Track>, double> >) simToGlbMuColl[simRef];
      rt = simToGlbMuColl[simRef];
      if (rt.size()!=0) {
	LogTrace(theCategory) << "\nTrackingParticle #" << i 
			      << " with pt=" << sqrt(simRef->momentum().perp2()) 
			      << " associated with quality:" << rt.begin()->second <<"\n";
      }
    }else{
      LogTrace(theCategory) 
	<< "TrackingParticle #" << i
	<< " with pt,eta,phi: " 
	<< sqrt(simRef->momentum().perp2()) << " , "
	<< simRef->momentum().eta() << " , "
	<< simRef->momentum().phi() << " , "
	<< " NOT associated to any reco::Track" << "\n";
    }

    if (rt.size()==0) continue;

    reco::TrackRef glbMuTrackRef = rt.begin()->first.castTo<TrackRef >();
    //edm::RefToBase<reco::Track> recoRef = rt.begin()->first; 

    tpBitMap_.insert(simRef,thisBit);
    recoBitMap_.insert( glbMuTrackRef, thisBit );
    recoTpMap_.insert( glbMuTrackRef, simRef );
    //unsigned int thisBit2;   
    //testMap2.insert(std::make_pair(thisBit,thisBit2) );   
    
  }
  ///////////////////////  


  LogDebug(theCategory)<<"MuonColl size " << muonColl.size();  
  // Analyzer reco::Muon
  for(View<Muon>::const_iterator iMuon = muonColl.begin();
      iMuon != muonColl.end(); ++iMuon) {

    int muClass = 0;

    LogTrace(theCategory)<<"Looking at a new muon.....";
    if ( iMuon->isGlobalMuon() ) {      
      LogTrace(theCategory)<<"which isGlobalMuon";

      //TrackingParticleRef & theTrp ;

      const TrackRef glbTrack = iMuon->combinedMuon();
      const RefToBase<Track> glbTrackRB(glbTrack);

      unsigned int theBit = recoBitMap_[glbTrack];
      TrackingParticleRef muTp = recoTpMap_[glbTrack];
      
      LogTrace(theCategory)<<"that is available " <<  glbTrack.isAvailable() << " with bit " << theBit;
      
      classifier_.evaluate( glbTrackRB );
      
      LogTrace(theCategory) <<"      with flags "; 
      for (std::size_t index = 0; index < classifier_.flags().size(); ++index){
	if (classifier_.flags()[index])
	  LogTrace(theCategory) <<"                 " 
				<<  TrackCategories::Names[index]; 
      }

      //LogTrace(theCategory) << "  and kink " << kink(glbTrack);

      std::vector<std::pair<TrackingParticleRef,double> > tpRefV;
      if ( glbTrack.isAvailable() && glbMuToSimColl.find(glbTrackRB) != glbMuToSimColl.end() ) {//get TP
	tpRefV = glbMuToSimColl[glbTrackRB];

	LogTrace(theCategory)<<"Found tpRefV of size " << tpRefV.size();

	///
	if (theBit == 0) {
	  LogTrace(theCategory) << "Sanity Check . . . . .";
	  std::vector<std::pair<TrackingParticleRef,double> >::const_iterator iIter;
	  std::vector<std::pair<TrackingParticleRef,double> >::const_iterator iBegin = tpRefV.begin();
	  std::vector<std::pair<TrackingParticleRef,double> >::const_iterator iEnd = tpRefV.end();
	  for( iIter = iBegin; iIter != iEnd; ++iIter) {
	    unsigned int iBit = getBit(iIter->first);
	    LogTrace(theCategory) << " . . . " << iBit;
	  }
	}
	///

	const TrackingParticleRef & trp = tpRefV.begin()->first;

	if (isPrimaryMuon(theBit)) muClass = 1;
	else if (isSiliconMuon(theBit)) muClass = 2;
	else if (isCalConversionMuon(theBit)) muClass = 3;
	else if (isOtherMuon(theBit)) muClass = 4;
	else muClass = 5;


	MuonME * thisME = 0;
	
	if(muClass==1) thisME = aME_;
	if(muClass==2) thisME = bME_;
	if(muClass==3) thisME = cME_;
	if(muClass==4) thisME = dME_;
	if(muClass==5) thisME = eME_;

	GlobalPoint pos(0,0,0);
	GlobalPoint decayPos(1,1,1);
	
	thisME->fill(*iMuon,trp,muTp);
	

      }//get TP


      LogTrace(theCategory)<<"MuClass " << muClass;


    }// isGlobal()
  }// loop over muon
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
GlbSelectorStudy::beginJob(const edm::EventSetup& iSetup)
{
  edm::Service<TFileService> fs;

  iSetup.getData( pdt_ );
//TFileDirectory prompt = fs->mkdir( "prompt" );
//decay = fs->mkdir( "decay" );
//outside = fs->mkdir( "outside" );
 // TH1F * h_pt = outside.make<TH1F>( "pt"  , "p_{t}", 100,  0., 100. );

 aME_ = new MuonME;
 bME_ = new MuonME;
 cME_ = new MuonME;
 dME_ = new MuonME;
 eME_ = new MuonME;
 //fME_ = new MuonME;
 //gME_ = new MuonME;
 // decayME_ = new MuonME;
 // outsideME_ = new MuonME;
 //
 aME_->bookHistograms(fs,"PrimaryMu");
 bME_->bookHistograms(fs,"SiliconMu");
 cME_->bookHistograms(fs,"CalConversionMu");
 dME_->bookHistograms(fs,"OtherMu");
 eME_->bookHistograms(fs,"PunchThrough");
 //fME_->bookHistograms(fs,"F");
 //gME_->bookHistograms(fs,"G");



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
  LogDebug(theCategory);
  typedef std::vector<Trajectory> TC;
  LogDebug(theCategory);
  TC staTrajs = theTrackTransformer->transform(muonIn);
  LogDebug(theCategory);
  if (staTrajs.empty()) {
    LogTrace(theCategory) << "Transformer: Add Traj failed!";
    return 0.;
  } 
  Trajectory muon(staTrajs.front());
  LogDebug(theCategory);
  
  double result = 0.0;
  
  typedef TransientTrackingRecHit::ConstRecHitPointer 	ConstRecHitPointer;
  typedef ConstRecHitPointer RecHit;
  typedef vector<TrajectoryMeasurement>::const_iterator TMI;
  LogDebug(theCategory);
  vector<TrajectoryMeasurement> meas = muon.measurements();
  LogDebug(theCategory);
  for ( TMI m = meas.begin(); m != meas.end(); m++ ) {
    LogDebug(theCategory);
    RecHit rhit = (*m).recHit();
    bool ok = false;
    if ( rhit->isValid() ) {
      //const GeomDetType& type = rhit->det()->detUnits().front()->type();
      //if ( type.module() == pixel || type.module() == silicon ) ok = true;
      
      if(DetId::Tracker == rhit->geographicalId().det()) ok = true;
      LogDebug(theCategory);
    }
    if ( !ok ) continue;
    
    const TrajectoryStateOnSurface& tsos = (*m).predictedState();
    LogDebug(theCategory);
    if ( tsos.isValid() ) {
      LogDebug(theCategory);
      double phi1 = tsos.globalPosition().phi();
      if ( phi1 < 0 ) phi1 = 2*M_PI + phi1;
      
      double phi2 = rhit->globalPosition().phi();
      if ( phi2 < 0 ) phi2 = 2*M_PI + phi2;
      
      double diff = fabs(phi1 - phi2);
      if ( diff > M_PI ) diff = 2*M_PI - diff;
      LogDebug(theCategory);
      //GlobalPoint hitPos = rhit->det()->toGlobal(rhit->localPosition());
      //GlobalError hitErr = rhit->det()->toGlobal(rhit->localPositionError());
      GlobalPoint hitPos = rhit->globalPosition();
      GlobalError hitErr = rhit->globalPositionError();
      LogDebug(theCategory);
      double error = hitErr.phierr(hitPos);  // error squared
      
      double s = ( error > 0.0 ) ? (diff*diff)/error : (diff*diff);
      result += s;
      LogDebug(theCategory);
    }
  }
  LogDebug(theCategory);
  return result;
  
}

void GlbSelectorStudy::printTruth(TrackingParticleRef & simRef, int i) {
  const TrackingParticle* simTP = simRef.get();
  TrackHistory const & tracer = classifier_.history();
  
  classifier_.evaluate( simRef );
  
    LogTrace(theCategory)<<"TP " << i << " with pdgID " << simTP->pdgId();
    LogTrace(theCategory)<<"     with GenSize " << simTP->genParticle().size();
    LogTrace(theCategory)<<"     with GeantSize " << simTP->g4Tracks().size();
    LogTrace(theCategory)<<"     with Pt      " << simTP->pt();
    LogTrace(theCategory)<<"     with Eta     " << simTP->eta();
    LogTrace(theCategory)<<"     with Phi     " << simTP->phi();
    LogTrace(theCategory)<<"     with Vtx Pos " << simTP->parentVertex()->position()
			 <<"       at " << sqrt(simTP->vertex().perp2()) << " , " <<  simTP->vertex().z();

    for(tv_iterator dv = simTP->decayVertices_begin(); dv != simTP->decayVertices_end(); ++dv) {
      LogTrace(theCategory)<<"      DecayVertex "<<(*dv)->position();
    }

    LogTrace(theCategory) <<"     with flags "; 
    for (std::size_t index = 0; index < classifier_.flags().size(); ++index){
      if (classifier_.flags()[index])
	LogTrace(theCategory) <<"                 " 
			      <<  TrackCategories::Names[index]; 
    }
    
    // Get the list of TrackingParticles associated to
    TrackHistory::SimParticleTrail simParticles(tracer.simParticleTrail());
    
    // Loop over all simParticles
    for (std::size_t hindex=0; hindex<simParticles.size(); hindex++)
      {
	LogTrace(theCategory)
	  << "  simParticles [" << hindex << "] : "
	  << particleString(simParticles[hindex]->pdgId())
	  //<< simParticles[hindex]->pdgId()
	  << " pT " <<  simParticles[hindex]->pt()
	  << " eta " <<  simParticles[hindex]->eta()
	  << " phi " <<  simParticles[hindex]->phi()
	  << " at " << sqrt(simParticles[hindex]->vertex().perp2()) << " , " <<  simParticles[hindex]->vertex().z() 
	  ;
      }
    
    // Get the list of TrackingVertexes associated to
    TrackHistory::SimVertexTrail simVertexes(tracer.simVertexTrail());
    
    // Loop over all simVertexes
    if ( !simVertexes.empty() )
      {
	for (std::size_t hindex=0; hindex<simVertexes.size(); hindex++)
	  {
	    LogTrace(theCategory) 
	      << "  simVertex    [" << hindex << "] : "
	      << vertexString(
			      simVertexes[hindex]->sourceTracks(),
			      simVertexes[hindex]->daughterTracks()
			      )
	      << " at " << simVertexes[hindex]->position().Rho() << " , " <<  simVertexes[hindex]->position().z() << " " <<  simVertexes[hindex]->inVolume()
	      ;
	  }
      }
    else
      LogTrace(theCategory) << "  simVertex not found" ;
    
    //    if(abs(simTP->pdgId()) == 13)
    //      for(TrackingParticle::g4t_iterator isimtk = simTP->g4Track_begin();isimtk!=simTP->g4Track_end();isimtk++) {//loop over sim in TP
    //	LogTrace(theCategory)<<"... look for SimMu MOTHER ....";
    //	MotherSearch mother(&*isimtk, simTracks, simVertexs, hepmc);
    //      }
  
}

std::string GlbSelectorStudy::particleString(int pdgId) const
{
    ParticleData const * pid;

    std::ostringstream vDescription;

    HepPDT::ParticleID particleType(pdgId);

    if (particleType.isValid())
    {
        pid = pdt_->particle(particleType);
        if (pid)
            vDescription << pid->name();
        else
            vDescription << pdgId;
    }
    else
        vDescription << pdgId;

    return vDescription.str();
}


std::string GlbSelectorStudy::vertexString(
    TrackingParticleRefVector in,
    TrackingParticleRefVector out
) const
{
    ParticleData const * pid;

    std::ostringstream vDescription;

    for (std::size_t j = 0; j < in.size(); j++)
    {
        if (!j) vDescription << "(";

        HepPDT::ParticleID particleType(in[j]->pdgId());

        if (particleType.isValid())
        {
            pid = pdt_->particle(particleType);
            if (pid)
                vDescription << pid->name();
            else
                vDescription << in[j]->pdgId();
        }
        else
            vDescription << in[j]->pdgId();

        if (j == in.size() - 1) vDescription << ")";
        else vDescription << ",";
    }

    vDescription << "->";

    for (std::size_t j = 0; j < out.size(); j++)
    {
        if (!j) vDescription << "(";

        HepPDT::ParticleID particleType(out[j]->pdgId());

        if (particleType.isValid())
        {
            pid = pdt_->particle(particleType);
            if (pid)
                vDescription << pid->name();
            else
                vDescription << out[j]->pdgId();
        }
        else
            vDescription << out[j]->pdgId();

        if (j == out.size() - 1) vDescription << ")";
        else vDescription << ",";
    }

    return vDescription.str();
}

unsigned int GlbSelectorStudy::getBit(const TrackingParticleRef& simRef) const{
    unsigned int thisBit = 0;
    
    if ( tpSelector_primary(*simRef) ) thisBit = primaryMuon;
    if ( tpSelector_silicon(*simRef) && ! isPrimaryMuon(thisBit) ) thisBit = siliconMuon;
    if ( tpSelector_calConversion(*simRef) && ! isSiliconMuon(thisBit) && ! isPrimaryMuon(thisBit) ) thisBit = calConversionMuon;
    if ( fabs(simRef->pdgId())==13 && ! isPrimaryMuon(thisBit) && ! isSiliconMuon(thisBit) && ! isCalConversionMuon(thisBit)) thisBit = otherMuon;

    return thisBit;
}

//define this as a plug-in
DEFINE_FWK_MODULE(GlbSelectorStudy);