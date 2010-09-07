#include "UserCode/Examples/plugins/PatTupleDumper.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "MuonAnalysis/Examples/interface/muonStations.h"

// for "luminosity"
#include "DataFormats/Common/interface/MergeableCounter.h"

// for selection cut
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TObjString.h>
#include <TDirectory.h>

#include <map>
#include <string>

#include "boost/lexical_cast.hpp"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;
using namespace edm;
using namespace reco;

using boost::lexical_cast;

/// Constructor
PatTupleDumper::PatTupleDumper(const edm::ParameterSet& pset):
    muons_(pset.getParameter<edm::InputTag>("muons")),
    selector_(pset.getParameter<std::string>("selection")),
    primaryVertices_(pset.getParameter<edm::InputTag>("primaryVertices")),
    outputFileName_(pset.getParameter<std::string>("outputFileName"))
{
  //edm::Service<TFileService> fs;

  //TFileDirectory md = fs->mkdir("metadata");
  //TDirectory *md_dir = md.cd();
  //md_dir->WriteTObject(new TObjString(muons_.encode().c_str()), "muons");
  //md_dir->WriteTObject(new TObjString(pset.getParameter<std::string>("selection").c_str()), "selection");
  ntupleInit();

}

/// Destructor
PatTupleDumper::~PatTupleDumper()
{
  delete outputFile_;
}

void PatTupleDumper::beginJob( ) {
  outputFile_  = new TFile( outputFileName_.c_str(), "RECREATE" ) ;
  
  myTree = new TTree("myTree", "myTree");

  myTree->Branch("evtRun",&evtRun,"evtRun/I");
  myTree->Branch("evtNum",&evtNum,"evtNum/I");  
  myTree->Branch("evtWeight",&evtWeight,"evtWeight/F");
  myTree->Branch("evtProcessId",&evtProcessId,"evtProcessId/I");
  
  
  myTree->Branch( "mSize",&mSize,"mSize/I");
  myTree->Branch( "mPass",&mPass,"mPass/I");
  myTree->Branch( "p",&p);
  myTree->Branch( "pt",&pt); 
  myTree->Branch( "eta",&eta); 
  myTree->Branch( "phi",&phi); 
  myTree->Branch( "charge",&charge); 
  
  myTree->Branch( "pSta",   &pSta); 
  myTree->Branch( "ptSta",  &ptSta); 
  myTree->Branch( "etaSta", &etaSta); 
  myTree->Branch( "phiSta", &phiSta); 
  
  myTree->Branch( "dxy",&dxy);
  myTree->Branch( "dz",&dz);
  
  myTree->Branch( "pixelHits",&pixelHits);
  myTree->Branch( "pixelLayers",&pixelLayers);
  myTree->Branch( "trackerHits",&trackerHits);
  myTree->Branch( "trackerLostHitsInner",  &trackerLostHitsInner);
  myTree->Branch( "trackerLostHitsMiddle", &trackerLostHitsMiddle);
  myTree->Branch( "trackerLostHitsOuter",  &trackerLostHitsOuter);
  myTree->Branch( "muonHits",&muonHits);
  myTree->Branch( "muonBadHits",&muonBadHits);
  myTree->Branch( "globalHits",&globalHits);
  myTree->Branch( "globalMuonHits",&globalMuonHits);
  myTree->Branch( "trackerChi2n",&trackerChi2n);
  myTree->Branch( "muonChi2n",&muonChi2n);
  myTree->Branch( "trackerChi2Rel",&trackerChi2Rel);
  myTree->Branch( "muonChi2Rel",&muonChi2Rel);
  myTree->Branch( "globalChi2n",&globalChi2n);
  
  myTree->Branch( "deltaPt",&deltaPt);
  myTree->Branch( "deltaPtn",&deltaPtn);
  
  myTree->Branch( "muonHitCounts", &muonHitCounts);
  myTree->Branch( "muonHitCountsAny", &muonHitCountsAny);

  /*
    for(size_t j =0; j<4; ++j) {
    std::string intLabel = lexical_cast<std::string>(j+1);
    myTree->Branch( "muonHitCounts"+intLabel+"any", "muonStationHits");
    myTree->Branch( "muonHitCountsv"+intLabel, "muonStationHits");
    myTree->Branch( "muonHitCountsdt"+intLabel+"any", "muonStationHits");
    myTree->Branch( "muonHitCountscsc"+intLabel+"any", "muonStationHits");
    myTree->Branch( "muonHitCountsrpc"+intLabel+"any", "muonStationHits");
    myTree->Branch( "muonHitCountsdt"+intLabel, "muonStationHits");
    myTree->Branch( "muonHitCountscsc"+intLabel, "muonStationHits");
    myTree->Branch( "muonHitCountsrpc"+intLabel, "muonStationHits");
    } 
  */

  myTree->Branch("muonHitCounts_1_any",&muonHitCounts_1_any);
  myTree->Branch("muonHitCounts_2_any",&muonHitCounts_2_any);
  myTree->Branch("muonHitCounts_3_any",&muonHitCounts_3_any);
  myTree->Branch("muonHitCounts_4_any",&muonHitCounts_4_any);

  myTree->Branch("muonHitCounts_1_valid",&muonHitCounts_1_valid);
  myTree->Branch("muonHitCounts_2_valid",&muonHitCounts_2_valid);
  myTree->Branch("muonHitCounts_3_valid",&muonHitCounts_3_valid);
  myTree->Branch("muonHitCounts_4_valid",&muonHitCounts_4_valid);

  myTree->Branch( "muonHitCountsratio",&muonHitCountsratio);
  myTree->Branch( "muonHitCountsrpcratio",&muonHitCountsrpcratio);
  
  myTree->Branch( "trackIso05", &trackIso05);
  myTree->Branch( "ecalIso05",  &ecalIso05);
  myTree->Branch( "hcalIso05",  &hcalIso05);
  myTree->Branch( "trackIso03", &trackIso03);
  myTree->Branch( "ecalIso03",  &ecalIso03);
  myTree->Branch( "hcalIso03",  &hcalIso03);
  myTree->Branch( "combRelIso03", &combRelIso03);
  myTree->Branch( "combRelIso05", &combRelIso05);
  
  myTree->Branch( "muonStationsValid",&muonStationsValid);
  myTree->Branch( "muonStationsAny",&muonStationsAny);
  myTree->Branch( "muonStationsDTValid",&muonStationsDTValid);
  myTree->Branch( "muonStationsDTAny",&muonStationsDTAny);
  myTree->Branch( "muonStationsCSCValid",&muonStationsCSCValid);
  myTree->Branch( "muonStationsCSCAny",&muonStationsCSCAny);
  myTree->Branch( "muonStationsRPCValid",&muonStationsRPCValid);
  myTree->Branch( "muonStationsRPCAny",&muonStationsRPCAny);
  myTree->Branch( "numberOfChambers",&numberOfChambers);
  myTree->Branch( "segmentMatchesArb_MaxDepth",&segmentMatchesArb_MaxDepth); 
  myTree->Branch( "segmentMatchesArb",&segmentMatchesArb); 
  myTree->Branch( "segmentMatchesArb_1",&segmentMatchesArb_1); 
  myTree->Branch( "segmentMatchesArb_2",&segmentMatchesArb_2); 
  myTree->Branch( "segmentMatchesArb_3",&segmentMatchesArb_3); 
  myTree->Branch( "segmentMatchesArb_4",&segmentMatchesArb_4); 
  myTree->Branch( "segmentMatchesNoArb",&segmentMatchesNoArb); 
  myTree->Branch( "segmentMatchesNoArb_1",&segmentMatchesNoArb_1); 
  myTree->Branch( "segmentMatchesNoArb_2",&segmentMatchesNoArb_2); 
  myTree->Branch( "segmentMatchesNoArb_3",&segmentMatchesNoArb_3); 
  myTree->Branch( "segmentMatchesNoArb_4",&segmentMatchesNoArb_4); 
  myTree->Branch( "segmentMatchesFailArb",&segmentMatchesFailArb); 
  myTree->Branch( "segmentCompatArb",&segmentCompatArb); 
  myTree->Branch( "segmentCompatNoArb",&segmentCompatNoArb); 
  myTree->Branch( "caloCompat",&caloCompat); 
  
  myTree->Branch( "TMOneStationAngLoose",&TMOneStationAngLoose);
  myTree->Branch( "TMLastStationLoose",&TMLastStationLoose);
  myTree->Branch( "AllGlobalMuons",&AllGlobalMuons);
  myTree->Branch( "AllTrackerMuons",&AllTrackerMuons);
  myTree->Branch( "GlobalMuonPromptTight",&GlobalMuonPromptTight);

  
  //myTree->Branch( "prodz", "z"); 
  //myTree->Branch( "prodr", "r");
  //myTree->Branch( "prodd", "r");
  //book2d(*fs,pset,"prodrz","rz");

  /*  
  if (pset.existsAs<edm::InputTag>("normalization")) {
    normalization_ = pset.getParameter<edm::InputTag>("normalization");
    luminosity = fs->make<TH1D>("normalization", "normalization", 1, 0, 1);
    luminosity->Sumw2();
  }
  */

}

void PatTupleDumper::endJob() {

	outputFile_->Write() ;
	outputFile_->Close() ;

}

//void PatTupleDumper::book(const TTree &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) 
//{
//
//}

void PatTupleDumper::ntupleInit(){
  p	 = new vector<float> ;
  pt       = new vector<float> ;
  eta       = new vector<float> ;
  phi       = new vector<float> ;
  charge       = new vector<float> ;
  pSta       = new vector<float> ;
  ptSta       = new vector<float> ;
  etaSta       = new vector<float> ;
  phiSta       = new vector<float> ;
  dxy       = new vector<float> ;
  dz       = new vector<float> ;
  pixelHits       = new vector<int> ;
  pixelLayers       = new vector<int> ;
  trackerHits       = new vector<int> ;
  trackerLostHitsInner     = new vector<int> ;
  trackerLostHitsMiddle      = new vector<int> ;
  trackerLostHitsOuter       = new vector<int> ;
  muonHits       = new vector<int> ;
  muonBadHits       = new vector<int> ;
  globalHits       = new vector<int> ;
  globalMuonHits       = new vector<int> ;
  trackerChi2n       = new vector<float> ;
  muonChi2n       = new vector<float> ;
  trackerChi2Rel       = new vector<float> ;
  muonChi2Rel       = new vector<float> ;
  globalChi2n       = new vector<float> ;
  deltaPt       = new vector<float> ;
  deltaPtn       = new vector<float> ;
  muonHitCounts       = new vector<int> ;
  muonHitCountsAny       = new vector<int> ;

  muonHitCounts_1_any = new vector<int>;
  muonHitCounts_2_any = new vector<int>;
  muonHitCounts_3_any = new vector<int>;
  muonHitCounts_4_any = new vector<int>;

  muonHitCounts_1_valid = new vector<int>;
  muonHitCounts_2_valid = new vector<int>;
  muonHitCounts_3_valid = new vector<int>;
  muonHitCounts_4_valid = new vector<int>;

  muonHitCountsratio       = new vector<float> ;
  muonHitCountsrpcratio       = new vector<float> ;
  trackIso05       = new vector<float> ;
  ecalIso05       = new vector<float> ;
  hcalIso05       = new vector<float> ;
  trackIso03       = new vector<float> ;
  ecalIso03       = new vector<float> ;
  hcalIso03       = new vector<float> ;
  combRelIso03       = new vector<float> ;
  combRelIso05       = new vector<float> ;
  muonStationsValid       = new vector<int> ;
  muonStationsAny       = new vector<int> ;
  muonStationsDTValid       = new vector<int> ;
  muonStationsDTAny       = new vector<int> ;
  muonStationsCSCValid       = new vector<int> ;
  muonStationsCSCAny       = new vector<int> ;
  muonStationsRPCValid       = new vector<int> ;
  muonStationsRPCAny       = new vector<int> ;
  numberOfChambers       = new vector<int> ;
  segmentMatchesArb_MaxDepth       = new vector<int> ;
  segmentMatchesArb       = new vector<int> ;
  segmentMatchesArb_1       = new vector<int> ;
  segmentMatchesArb_2       = new vector<int> ;
  segmentMatchesArb_3       = new vector<int> ;
  segmentMatchesArb_4       = new vector<int> ;
  segmentMatchesNoArb       = new vector<int> ;
  segmentMatchesNoArb_1       = new vector<int> ;
  segmentMatchesNoArb_2       = new vector<int> ;
  segmentMatchesNoArb_3       = new vector<int> ;
  segmentMatchesNoArb_4       = new vector<int> ;
  segmentMatchesFailArb       = new vector<int> ;
  segmentCompatArb       = new vector<float> ;
  segmentCompatNoArb       = new vector<float> ;
  caloCompat       = new vector<float> ;
  TMOneStationAngLoose       = new vector<int> ;
  TMLastStationLoose       = new vector<int> ;
  AllGlobalMuons       = new vector<int> ;
  AllTrackerMuons       = new vector<int> ;
  GlobalMuonPromptTight       = new vector<int> ;

}

void PatTupleDumper::analyze(const edm::Event & event, const edm::EventSetup& eventSetup){
    using namespace edm;
    using namespace std;

    evtRun = (int)event.id().run();
    evtNum = (int)event.id().event();

    Handle<View<reco::Muon> > muons;
    event.getByLabel(muons_, muons);
    mSize = (int)muons->size();
    LogDebug("Dumper")<<"Muons size: "<<muons->size();

    Handle<vector<reco::Vertex> > vertices;
    event.getByLabel(primaryVertices_, vertices);
    int k =0;
    mPass=0;
    foreach (const reco::Muon &recomu, *muons) {
      // we want to make a pat::Muon so that we can access directly muonID in the cuts
      
      const pat::Muon &mu = (typeid(recomu) == typeid(pat::Muon) ? static_cast<const pat::Muon &>(recomu) : pat::Muon(recomu));


      if (!selector_(mu)) continue;
      mPass++;
      p->push_back(mu.p());
      pt->push_back(mu.pt());
      eta->push_back(mu.eta());
      phi->push_back(mu.phi());
      charge->push_back(mu.charge());
      
      TMOneStationAngLoose->push_back(muon::isGoodMuon(mu,muon::TMOneStationAngLoose));
      TMLastStationLoose->push_back(muon::isGoodMuon(mu,muon::TMLastStationLoose));
      AllGlobalMuons->push_back(muon::isGoodMuon(mu,muon::AllGlobalMuons));
      AllTrackerMuons->push_back(muon::isGoodMuon(mu,muon::AllTrackerMuons));
      GlobalMuonPromptTight->push_back(muon::isGoodMuon(mu,muon::GlobalMuonPromptTight));
      
      if (mu.innerTrack().isNonnull()) {
	pixelHits->push_back(mu.innerTrack()->hitPattern().numberOfValidPixelHits());
	pixelLayers->push_back(mu.innerTrack()->hitPattern().pixelLayersWithMeasurement());
	trackerHits->push_back(mu.innerTrack()->hitPattern().numberOfValidHits()); 
	trackerLostHitsMiddle->push_back(mu.innerTrack()->hitPattern().numberOfLostHits());
	trackerLostHitsInner->push_back(mu.innerTrack()->trackerExpectedHitsInner().numberOfLostHits());
	trackerLostHitsOuter->push_back(mu.innerTrack()->trackerExpectedHitsOuter().numberOfLostHits());
	trackerChi2n->push_back(mu.innerTrack()->normalizedChi2());
	
	if (!vertices->empty() && !vertices->front().isFake()) {
	  const reco::Vertex &vtx = vertices->front();
	  dxy->push_back(mu.innerTrack()->dxy(vtx.position()));
	  dz->push_back(mu.innerTrack()->dz(vtx.position()));
	}
      }
      
      if (mu.outerTrack().isNonnull()) {
	pSta->push_back(mu.outerTrack()->p());
	ptSta->push_back(mu.outerTrack()->pt());
	etaSta->push_back(mu.outerTrack()->eta());
	phiSta->push_back(mu.outerTrack()->phi());
	
	muonHits->push_back(mu.outerTrack()->numberOfValidHits());
	muonBadHits->push_back(mu.outerTrack()->recHitsSize() - mu.outerTrack()->numberOfValidHits());
	muonChi2n->push_back(mu.outerTrack()->normalizedChi2());
	
	if(mu.hasUserInt("muonStations")) {
	  muonStationsAny->push_back(mu.userInt("muonStations:any"));
	  muonStationsValid->push_back(mu.userInt("muonStations"));
	  muonStationsDTAny->push_back(mu.userInt("muonStations:dtAny"));
	  muonStationsDTValid->push_back(mu.userInt("muonStations:dt"));
	  muonStationsCSCAny->push_back(mu.userInt("muonStations:cscAny"));
	  muonStationsCSCValid->push_back(mu.userInt("muonStations:csc"));
	  muonStationsRPCAny->push_back(mu.userInt("muonStations:rpcAny"));
	  muonStationsRPCValid->push_back(mu.userInt("muonStations:rpc"));
	} else if ( ( mu.outerTrack()->extra().isAvailable()   ) && 
		    ( mu.outerTrack()->recHitsSize() > 0       ) &&
		    ( mu.outerTrack()->recHit(0).isAvailable() )     ) {
	  
	  muonStationsValid->push_back(muon::muonStations(mu.outerTrack(), 0, true));
	  muonStationsAny->push_back(muon::muonStations(mu.outerTrack(), 0, false));
	  //if(mu.hasUserInt("muonStations:Any"));
	  //muonStationsAny->push_back(mu.userInt("muonStations:Any"));
	  float abseta = std::abs(mu.outerTrack()->eta());
	  if (abseta <= 1.2) {
	    muonStationsDTValid->push_back(muon::muonStations(mu.outerTrack(),MuonSubdetId::DT, true));
	    muonStationsDTAny->push_back(muon::muonStations(mu.outerTrack(),MuonSubdetId::DT, false));
	  } 
	  if (abseta <= 1.6) {
	    muonStationsRPCValid->push_back(muon::muonStations(mu.outerTrack(),MuonSubdetId::RPC, true));
	    muonStationsRPCAny->push_back(muon::muonStations(mu.outerTrack(),MuonSubdetId::RPC, false));
	  } 
	  if (abseta >= 0.8) {
	    muonStationsCSCValid->push_back(muon::muonStations(mu.outerTrack(),MuonSubdetId::CSC, true));
	    muonStationsCSCAny->push_back(muon::muonStations(mu.outerTrack(),MuonSubdetId::CSC, false));
	  }
	}
      }
      
      if (mu.globalTrack().isNonnull()) {
	globalHits->push_back(mu.globalTrack()->numberOfValidHits());
	globalMuonHits->push_back(mu.globalTrack()->hitPattern().numberOfValidMuonHits());
	globalChi2n->push_back(mu.globalTrack()->normalizedChi2());
	deltaPt->push_back((mu.outerTrack()->pt() - mu.innerTrack()->pt()));
	deltaPtn->push_back((mu.outerTrack()->pt() - mu.innerTrack()->pt())/(mu.innerTrack()->pt()));
      }
      
      if(mu.isQualityValid()) {
	muonChi2Rel->push_back(mu.combinedQuality().staRelChi2);
	trackerChi2Rel->push_back(mu.combinedQuality().trkRelChi2);
      }
      
      if (mu.isIsolationValid()) {
	trackIso05->push_back(mu.isolationR05().sumPt);
	ecalIso05->push_back(mu.isolationR05().emEt);
	hcalIso05->push_back(mu.isolationR05().hadEt);
	trackIso03->push_back(mu.isolationR03().sumPt);
	ecalIso03->push_back(mu.isolationR03().emEt);
	hcalIso03->push_back(mu.isolationR03().hadEt);
	combRelIso03->push_back( (mu.isolationR03().sumPt + mu.isolationR03().emEt + mu.isolationR03().hadEt) / mu.pt() );
	combRelIso05->push_back( (mu.isolationR05().sumPt + mu.isolationR05().emEt + mu.isolationR05().hadEt) / mu.pt() );
      }
      
      if (mu.isMatchesValid()) {
	
	numberOfChambers->push_back(mu.numberOfChambers());
	
	segmentMatchesArb->push_back(mu.numberOfMatches(reco::Muon::SegmentAndTrackArbitration));
	segmentMatchesNoArb->push_back(mu.numberOfMatches(reco::Muon::SegmentArbitration));
	segmentMatchesFailArb->push_back(mu.numberOfMatches(reco::Muon::SegmentArbitration) - mu.numberOfMatches(reco::Muon::SegmentAndTrackArbitration));
	
	//adam stations with matched segments
	unsigned int maskST_Arb = mu.stationMask(reco::Muon::SegmentAndTrackArbitration);
	unsigned int maskS_Arb = mu.stationMask(reco::Muon::SegmentArbitration);
	
	int maxDepth = 0;
	if((maskST_Arb & 1<<0)||(maskST_Arb & 1<<4)) maxDepth = 1;
	if((maskST_Arb & 1<<1)||(maskST_Arb & 1<<5)) maxDepth = 2;
	if((maskST_Arb & 1<<2)||(maskST_Arb & 1<<6)) maxDepth = 3;
	if((maskST_Arb & 1<<3)||(maskST_Arb & 1<<7)) maxDepth = 4;
	
	segmentMatchesArb_MaxDepth->push_back(maxDepth);
	
	segmentMatchesArb_1->push_back(((maskST_Arb & 1<<0)||(maskST_Arb & 1<<4)));
	segmentMatchesArb_2->push_back(((maskST_Arb & 1<<1)||(maskST_Arb & 1<<5)));
	segmentMatchesArb_3->push_back(((maskST_Arb & 1<<2)||(maskST_Arb & 1<<6)));
	segmentMatchesArb_4->push_back(((maskST_Arb & 1<<3)||(maskST_Arb & 1<<7)));
	
	segmentMatchesNoArb_1->push_back(((maskS_Arb & 1<<0)||(maskS_Arb & 1<<4)));
	segmentMatchesNoArb_2->push_back(((maskS_Arb & 1<<1)||(maskS_Arb & 1<<5)));
	segmentMatchesNoArb_3->push_back(((maskS_Arb & 1<<2)||(maskS_Arb & 1<<6)));
	segmentMatchesNoArb_4->push_back(((maskS_Arb & 1<<3)||(maskS_Arb & 1<<7)));
	
	//adam end stations with matched segments
	
	segmentCompatArb->push_back(muon::segmentCompatibility(mu, reco::Muon::SegmentAndTrackArbitration));
	segmentCompatNoArb->push_back(muon::segmentCompatibility(mu, reco::Muon::SegmentArbitration));
      }
      
      if (mu.isCaloCompatibilityValid()) {
	caloCompat->push_back(mu.caloCompatibility());
      }
      
      if(mu.hasUserInt("muonHitCounts")) {
	
	muonHitCounts->push_back(mu.userInt("muonHitCounts"));
	muonHitCountsAny->push_back(mu.userInt("muonHitCounts:any"));
	
	std::string any = "any";
	std::string valid = "v";
	int rpcTotalHits = 0;
	
	muonHitCounts_1_any->push_back(mu.userInt("muonHitCounts:1"+any));
	muonHitCounts_2_any->push_back(mu.userInt("muonHitCounts:2"+any));
	muonHitCounts_3_any->push_back(mu.userInt("muonHitCounts:3"+any));
	muonHitCounts_4_any->push_back(mu.userInt("muonHitCounts:4"+any));
	
	muonHitCounts_1_valid->push_back(mu.userInt("muonHitCounts:v1"));
	muonHitCounts_2_valid->push_back(mu.userInt("muonHitCounts:v2"));
	muonHitCounts_3_valid->push_back(mu.userInt("muonHitCounts:v3"));
	muonHitCounts_4_valid->push_back(mu.userInt("muonHitCounts:v4"));
	
	
	for(size_t j =0; j<4; ++j) {
	  std::string intLabel = lexical_cast<std::string>(j+1);
	  rpcTotalHits += mu.userInt("muonHitCounts:rpc"+intLabel); 
	}
	
	//muonHitCounts_st_any->push_back(anyHit);
	/*
	  for(size_t j =0; j<4; ++j) {
	  std::string intLabel = lexical_cast<std::string>(j+1);
	  muonHitCounts"+intLabel+any]->Fill(mu.userInt("muonHitCounts:"+intLabel+any));
	  muonHitCounts"+valid+intLabel]->Fill(mu.userInt("muonHitCounts:"+valid+intLabel));
	  muonHitCountsdt"+intLabel+any]->Fill(mu.userInt("muonHitCounts:dt"+intLabel+any));
	  muonHitCountsdt"+intLabel]->Fill(mu.userInt("muonHitCounts:dt"+intLabel));
	  muonHitCountscsc"+intLabel+any]->Fill(mu.userInt("muonHitCounts:csc"+intLabel+any));
	  muonHitCountscsc"+intLabel]->Fill(mu.userInt("muonHitCounts:csc"+intLabel));
	  muonHitCountsrpc"+intLabel+any]->Fill(mu.userInt("muonHitCounts:rpc"+intLabel+any));
	  muonHitCountsrpc"+intLabel]->Fill(mu.userInt("muonHitCounts:rpc"+intLabel));
	  rpcTotalHits += mu.userInt("muonHitCounts:rpc"+intLabel);
	  } 
	*/
	
	if(mu.userInt("muonHitCounts") > 0) muonHitCountsratio->push_back(float(mu.userInt("muonHitCounts:v1"))/float(mu.userInt("muonHitCounts")));
	if(rpcTotalHits > 0) muonHitCountsrpcratio->push_back(float(mu.userInt("muonHitCounts:rpc1"))/float(rpcTotalHits));
	
      }
      
      float z = 0.;
      float rho = 0.;
      float d = 0.;
      
      if(mu.hasUserFloat("classByHitsGlb:prodZ")) z = mu.userFloat("classByHitsGlb:prodZ");
      if(mu.hasUserFloat("classByHitsGlb:prodRho")) rho = mu.userFloat("classByHitsGlb:prodRho");
      d = sqrt(z*z + rho*rho);
      
      //prodz->push_back(z);
      //prodr->push_back(rho);
      //prodd->push_back(d);
      //prodrz->push_back(z,rho);

      myTree->Fill();      
      
    }
    

    
}

/*
  void PatTupleDumper::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) 
  {
  if (luminosity != 0) {
  edm::Handle<edm::MergeableCounter> mc;
  iLumi.getByLabel(normalization_, mc);
  luminosity->Fill(0.5, double(mc->value));
  // set the correct uncertainty from counting statistics
  luminosity->SetBinError(1, sqrt(luminosity->GetBinContent(1)));
  }
  }
*/

DEFINE_FWK_MODULE(PatTupleDumper);







