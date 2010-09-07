#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"

//#include "FWCore/Framework/interface/TriggerNames.h"

class TFile;

#include <TMath.h>
#include "TFile.h"
#include <Math/VectorUtil.h>
#include <TVector3.h>
#include <utility>

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "MuonAnalysis/Examples/interface/muonStations.h"
//#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include <vector>
#include <TTree.h>

using namespace std;

class PatTupleDumper : public edm::EDAnalyzer {
 public:
  explicit PatTupleDumper( const edm::ParameterSet & );
  ~PatTupleDumper();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;
  //void book(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) ;
  //void book(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name) { book(fs,pset,name,name); }
  //float DeltaPhi(float v1, float v2);
  //float GetDeltaR(float eta1, float eta2, float phi1, float phi2);
  //void _fillGmr(const edm::Handle<reco::MuonCollection> muons);
  void ntupleInit();	
 private:

  edm::InputTag muons_;
  StringCutObjectSelector<pat::Muon> selector_;
  edm::InputTag primaryVertices_;
  std::string outputFileName_;

  edm::InputTag normalization_;
  TFile * outputFile_ ;
  TTree * myTree;
  // - - - -for Event
  int evtRun;
  int evtNum;
  float evtWeight;
  int evtProcessId;

  // - -  - -Deault Global Muon Collection
  int mSize ;
  int mPass ;
  vector<float> *p	      ;
  vector<float> *pt            ;
  vector<float> *eta            ;
  vector<float> *phi            ;
  vector<float> *charge            ;
  vector<float> *pSta            ;
  vector<float> *ptSta            ;
  vector<float> *etaSta            ;
  vector<float> *phiSta            ;
  vector<float> *dxy            ;
  vector<float> *dz            ;
  vector<int> *pixelHits            ;
  vector<int> *pixelLayers            ;
  vector<int> *trackerHits            ;
  vector<int> *trackerLostHitsInner          ;
  vector<int> *trackerLostHitsMiddle           ;
  vector<int> *trackerLostHitsOuter            ;
  vector<int> *muonHits            ;
  vector<int> *muonBadHits            ;
  vector<int> *globalHits            ;
  vector<int> *globalMuonHits            ;
  vector<float> *trackerChi2n            ;
  vector<float> *muonChi2n            ;
  vector<float> *trackerChi2Rel            ;
  vector<float> *muonChi2Rel            ;
  vector<float> *globalChi2n            ;
  vector<float> *deltaPt            ;
  vector<float> *deltaPtn            ;
  vector<int> *muonHitCounts            ;
  vector<int> *muonHitCountsAny            ;

  vector<int> *muonHitCounts_1_any;
  vector<int> *muonHitCounts_2_any;
  vector<int> *muonHitCounts_3_any;
  vector<int> *muonHitCounts_4_any;

  vector<int> *muonHitCounts_1_valid;
  vector<int> *muonHitCounts_2_valid;
  vector<int> *muonHitCounts_3_valid;
  vector<int> *muonHitCounts_4_valid;


  vector<float> *muonHitCountsratio            ;
  vector<float> *muonHitCountsrpcratio            ;
  vector<float> *trackIso05            ;
  vector<float> *ecalIso05            ;
  vector<float> *hcalIso05            ;
  vector<float> *trackIso03            ;
  vector<float> *ecalIso03            ;
  vector<float> *hcalIso03            ;
  vector<float> *combRelIso03            ;
  vector<float> *combRelIso05            ;
  vector<int> *muonStationsValid            ;
  vector<int> *muonStationsAny            ;
  vector<int> *muonStationsDTValid            ;
  vector<int> *muonStationsDTAny            ;
  vector<int> *muonStationsCSCValid            ;
  vector<int> *muonStationsCSCAny            ;
  vector<int> *muonStationsRPCValid            ;
  vector<int> *muonStationsRPCAny            ;
  vector<int> *numberOfChambers            ;
  vector<int> *segmentMatchesArb_MaxDepth            ;
  vector<int> *segmentMatchesArb            ;
  vector<int> *segmentMatchesArb_1            ;
  vector<int> *segmentMatchesArb_2            ;
  vector<int> *segmentMatchesArb_3            ;
  vector<int> *segmentMatchesArb_4            ;
  vector<int> *segmentMatchesNoArb            ;
  vector<int> *segmentMatchesNoArb_1            ;
  vector<int> *segmentMatchesNoArb_2            ;
  vector<int> *segmentMatchesNoArb_3            ;
  vector<int> *segmentMatchesNoArb_4            ;
  vector<int> *segmentMatchesFailArb            ;
  vector<float> *segmentCompatArb            ;
  vector<float> *segmentCompatNoArb            ;
  vector<float> *caloCompat            ;
  vector<int> *TMOneStationAngLoose            ;
  vector<int> *TMLastStationLoose            ;
  vector<int> *AllGlobalMuons            ;
  vector<int> *AllTrackerMuons            ;
  vector<int> *GlobalMuonPromptTight            ;
  //  vector<float> *            ;
  //  vector<float> *            ;
  //  vector<float> *            ;
  //  vector<float> *            ;
  //  vector<float> *            ;
  //  vector<float> *            ;

};
