/** \class InclusiveMuonPlotsMRTU
 *  Make inclusive muon plots
 *
 *  Based on r 1.4 of InclusiveMuonPlots by G.P.  This version 
 *  adds figures that only I am interested in . . . 
 *
 *  \author G. Petrucciani - UCSD (Giovanni.Petrucciani@cern.ch), ...
 */

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

using boost::lexical_cast;

class InclusiveMuonPlotsMRTU: public edm::EDAnalyzer {
    public:
        /// Constructor
        InclusiveMuonPlotsMRTU(const edm::ParameterSet& pset) ;

        /// Destructor
        virtual ~InclusiveMuonPlotsMRTU() ;

        // Operations
        void analyze(const edm::Event & event, const edm::EventSetup& eventSetup) ;

        void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup&);

        void book(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) ;
        void book2d(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) ;
        void book(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name) { book(fs,pset,name,name); }

        void bookProf(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) ;
        void bookProf(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name) { bookProf(fs,pset,name,name); }

    private:
        edm::InputTag muons_;
        StringCutObjectSelector<pat::Muon> selector_;
        StringCutObjectSelector<pat::Muon> subSelector_;

        bool old36Xdata_;

        edm::InputTag primaryVertices_;
        edm::InputTag normalization_;

        // we don't care too much about performance
        std::map<std::string, TH1*>      plots;
        std::map<std::string, TProfile*> profiles;

        TH1D *luminosity;
};

/// Constructor
InclusiveMuonPlotsMRTU::InclusiveMuonPlotsMRTU(const edm::ParameterSet& pset):
  muons_(pset.getParameter<edm::InputTag>("muons")),
  selector_(pset.getParameter<std::string>("selection")),
  subSelector_(pset.existsAs<std::string>("subSelection") ? pset.getParameter<std::string>("subSelection") : ""),
  old36Xdata_(pset.existsAs<bool>("old36Xdata") ? pset.getParameter<bool>("old36Xdata") : true),
  primaryVertices_(pset.getParameter<edm::InputTag>("primaryVertices")),
  luminosity(0) // by default, we don't have luminosity info
{
  edm::Service<TFileService> fs;
  
  TFileDirectory md = fs->mkdir("metadata");
  TDirectory *md_dir = md.cd();
  md_dir->WriteTObject(new TObjString(muons_.encode().c_str()), "muons");
  md_dir->WriteTObject(new TObjString(pset.getParameter<std::string>("selection").c_str()), "selection");
    
  book(*fs, pset, "trackerChi2Rel");
  book(*fs, pset, "muonChi2Rel");
  book(*fs, pset, "chi2LocalPosition");
  book(*fs, pset, "chi2LocalMomentum");
  book(*fs, pset, "localDistance");
  book(*fs, pset, "globalDeltaEtaPhi");
  book(*fs, pset, "tightMatch","bool");
  book(*fs, pset, "glbTrackProbability");
  book2d(*fs, pset, "chi2LocalPositionlocalDistance","chi2ld");
  book2d(*fs, pset, "chi2LocalMomentumlocalDistance","chi2mld");
  
  book(*fs, pset, "deltaPt");
  book(*fs, pset, "deltaPtn");
  
  book(*fs, pset, "muonHitCounts", "muonHits");
  book(*fs, pset, "muonHitCountsany", "muonHits");
  for(size_t j =0; j<4; ++j) {
    std::string intLabel = lexical_cast<std::string>(j+1);
    book(*fs, pset, "muonHitCounts"+intLabel+"any", "muonStationHits");
    book(*fs, pset, "muonHitCountsv"+intLabel, "muonStationHits");
    book(*fs, pset, "muonHitCountsdt"+intLabel+"any", "muonStationHits");
    book(*fs, pset, "muonHitCountscsc"+intLabel+"any", "muonStationHits");
    book(*fs, pset, "muonHitCountsrpc"+intLabel+"any", "muonStationHits");
    book(*fs, pset, "muonHitCountsdt"+intLabel, "muonStationHits");
    book(*fs, pset, "muonHitCountscsc"+intLabel, "muonStationHits");
    book(*fs, pset, "muonHitCountsrpc"+intLabel, "muonStationHits");
  }
  book(*fs, pset, "muonHitCountsratio","ratio");
  book(*fs, pset, "muonHitCountsrpcratio","ratio");
  
  book(*fs, pset, "numberOfChambers",     "segmentMatches");
  book(*fs, pset, "segmentMatchesArb_MaxDepth","segmentMatches"); 
  book(*fs, pset, "segmentMatchesArb_1",  "bool"); 
  book(*fs, pset, "segmentMatchesArb_2",  "bool"); 
  book(*fs, pset, "segmentMatchesArb_3",  "bool"); 
  book(*fs, pset, "segmentMatchesArb_4",  "bool"); 
  book(*fs, pset, "segmentMatchesNoArb_1","bool"); 
  book(*fs, pset, "segmentMatchesNoArb_2","bool"); 
  book(*fs, pset, "segmentMatchesNoArb_3","bool"); 
  book(*fs, pset, "segmentMatchesNoArb_4","bool"); 
  
  book(*fs, pset, "TMLastStationLoose", "bool");
  
  book(*fs, pset, "prodz", "z"); 
  book(*fs, pset, "prodr", "r");
  book(*fs, pset, "prodd", "r");
  book2d(*fs,pset,"prodrz","rz");
  
  if (pset.existsAs<edm::InputTag>("normalization")) {
    normalization_ = pset.getParameter<edm::InputTag>("normalization");
    luminosity = fs->make<TH1D>("normalization", "normalization", 1, 0, 1);
    luminosity->Sumw2();
  }
  
}

/// Destructor
InclusiveMuonPlotsMRTU::~InclusiveMuonPlotsMRTU()
{

}

void InclusiveMuonPlotsMRTU::book(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) 
{
    typedef std::vector<double> vdouble;
    if (pset.existsAs<vdouble>(basename+"Bins")) {
        vdouble bins = pset.getParameter<vdouble>(basename+"Bins");
        plots[name] = fs.make<TH1D>(name.c_str(), name.c_str(), bins.size()-1, &bins[0]);
    } else {
        uint32_t nbins = pset.getParameter<uint32_t>(basename+"Bins");
        vdouble  range = pset.getParameter<vdouble>(basename+"Range");
        if (range.size() != 2) throw cms::Exception("Configuration") << "parameter '" << basename << "Range' is not of the form (min, max).\n";
        plots[name] = fs.make<TH1D>(name.c_str(), name.c_str(), nbins, range[0], range[1]);
    }
}

void InclusiveMuonPlotsMRTU::book2d(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) 
{
    typedef std::vector<double> vdouble;
    //if (pset.existsAs<vdouble>(basename+"Bins")) {
    //vdouble bins = pset.getParameter<vdouble>(basename+"Bins");
    //plots[name] = fs.make<TH1D>(name.c_str(), name.c_str(), bins.size()-1, &bins[0]);
    //} else {
    uint32_t nbinsx = pset.getParameter<uint32_t>(basename+"XBins");
    vdouble  rangex = pset.getParameter<vdouble>(basename+"XRange");
    uint32_t nbinsy = pset.getParameter<uint32_t>(basename+"YBins");
    vdouble  rangey = pset.getParameter<vdouble>(basename+"YRange");
    if (rangex.size() != 2 || rangey.size() != 2) throw cms::Exception("Configuration") << "parameter '" << basename <<"Range' is not of the form (min, max).\n";
    plots[name] = fs.make<TH2D>(name.c_str(), name.c_str(), nbinsx, rangex[0], rangex[1],nbinsy,rangey[0],rangey[1]);
    //}
}

void InclusiveMuonPlotsMRTU::bookProf(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) 
{
    typedef std::vector<double> vdouble;
    if (pset.existsAs<vdouble>(basename+"Bins")) {
        vdouble bins = pset.getParameter<vdouble>(basename+"Bins");
        profiles[name] = fs.make<TProfile>(name.c_str(), name.c_str(), bins.size()-1, &bins[0]);
    } else {
        uint32_t nbins = pset.getParameter<uint32_t>(basename+"Bins");
        vdouble  range = pset.getParameter<vdouble>(basename+"Range");
        if (range.size() != 2) throw cms::Exception("Configuration") << "parameter '" << basename << "Range' is not of the form (min, max).\n";
        profiles[name] = fs.make<TProfile>(name.c_str(), name.c_str(), nbins, range[0], range[1]);
    }
}


void InclusiveMuonPlotsMRTU::analyze(const edm::Event & event, const edm::EventSetup& eventSetup){
    using namespace edm;
    using namespace std;

    Handle<View<reco::Muon> > muons;
    event.getByLabel(muons_, muons);

    Handle<vector<reco::Vertex> > vertices;
    event.getByLabel(primaryVertices_, vertices);

    size_t nmu = 0;
    foreach (const reco::Muon &recomu, *muons) {
      // we want to make a pat::Muon so that we can access directly muonID in the cuts
      const pat::Muon &mu = (typeid(recomu) == typeid(pat::Muon) ? static_cast<const pat::Muon &>(recomu) : pat::Muon(recomu));
      
      if (!selector_(mu)) continue;
      nmu++;
      if (!subSelector_(mu)) continue; // apply after counting
      
      plots["TMLastStationLoose"]->Fill(muon::isGoodMuon(mu,muon::TMLastStationLoose));
      
      if (mu.globalTrack().isNonnull()) {
	plots["deltaPt"]->Fill((mu.outerTrack()->pt() - mu.innerTrack()->pt()));
	plots["deltaPtn"]->Fill((mu.outerTrack()->pt() - mu.innerTrack()->pt())/(mu.innerTrack()->pt()));
      }
      
      if(mu.isQualityValid()) {
	plots["muonChi2Rel"]->Fill(mu.combinedQuality().staRelChi2);
	plots["trackerChi2Rel"]->Fill(mu.combinedQuality().trkRelChi2);
	plots["chi2LocalPosition"]->Fill(mu.combinedQuality().chi2LocalPosition);
	plots["chi2LocalMomentum"]->Fill(mu.combinedQuality().chi2LocalMomentum);
	plots["localDistance"]->Fill(mu.combinedQuality().localDistance);
	plots["globalDeltaEtaPhi"]->Fill(mu.combinedQuality().globalDeltaEtaPhi);
	plots["tightMatch"]->Fill(mu.combinedQuality().tightMatch);
	plots["glbTrackProbability"]->Fill(mu.combinedQuality().glbTrackProbability);
	plots["chi2LocalPositionlocalDistance"]->Fill(mu.combinedQuality().chi2LocalPosition,mu.combinedQuality().localDistance);
	plots["chi2LocalMomentumlocalDistance"]->Fill(mu.combinedQuality().chi2LocalMomentum,mu.combinedQuality().localDistance);
      }
      
      if (mu.isMatchesValid()) {	
	plots["numberOfChambers"]->Fill(mu.numberOfChambers());
	
	//adam stations with matched segments
	unsigned int maskST_Arb = mu.stationMask(reco::Muon::SegmentAndTrackArbitration);
	unsigned int maskS_Arb = mu.stationMask(reco::Muon::SegmentArbitration);
	
	int maxDepth = 0;
	if((maskST_Arb & 1<<0)||(maskST_Arb & 1<<4)) maxDepth = 1;
	if((maskST_Arb & 1<<1)||(maskST_Arb & 1<<5)) maxDepth = 2;
	if((maskST_Arb & 1<<2)||(maskST_Arb & 1<<6)) maxDepth = 3;
	if((maskST_Arb & 1<<3)||(maskST_Arb & 1<<7)) maxDepth = 4;
	
	plots["segmentMatchesArb_MaxDepth"]->Fill(maxDepth);
	
	plots["segmentMatchesArb_1"]->Fill(((maskST_Arb & 1<<0)||(maskST_Arb & 1<<4)));
	plots["segmentMatchesArb_2"]->Fill(((maskST_Arb & 1<<1)||(maskST_Arb & 1<<5)));
	plots["segmentMatchesArb_3"]->Fill(((maskST_Arb & 1<<2)||(maskST_Arb & 1<<6)));
	plots["segmentMatchesArb_4"]->Fill(((maskST_Arb & 1<<3)||(maskST_Arb & 1<<7)));
	
	plots["segmentMatchesNoArb_1"]->Fill(((maskS_Arb & 1<<0)||(maskS_Arb & 1<<4)));
	plots["segmentMatchesNoArb_2"]->Fill(((maskS_Arb & 1<<1)||(maskS_Arb & 1<<5)));
	plots["segmentMatchesNoArb_3"]->Fill(((maskS_Arb & 1<<2)||(maskS_Arb & 1<<6)));
	plots["segmentMatchesNoArb_4"]->Fill(((maskS_Arb & 1<<3)||(maskS_Arb & 1<<7)));
	
	//adam end stations with matched segments
	
        }

      
      if(mu.hasUserInt("muonHitCounts")) {
	
	plots["muonHitCounts"]->Fill(mu.userInt("muonHitCounts"));
	plots["muonHitCountsany"]->Fill(mu.userInt("muonHitCounts:any"));
	
	std::string any = "any";
	std::string valid = "v";
	int rpcTotalHits = 0;
	for(size_t j =0; j<4; ++j) {
	  std::string intLabel = lexical_cast<std::string>(j+1);
	  plots["muonHitCounts"+intLabel+any]->Fill(mu.userInt("muonHitCounts:"+intLabel+any));
	  plots["muonHitCounts"+valid+intLabel]->Fill(mu.userInt("muonHitCounts:"+valid+intLabel));
	  plots["muonHitCountsdt"+intLabel+any]->Fill(mu.userInt("muonHitCounts:dt"+intLabel+any));
	  plots["muonHitCountsdt"+intLabel]->Fill(mu.userInt("muonHitCounts:dt"+intLabel));
	  plots["muonHitCountscsc"+intLabel+any]->Fill(mu.userInt("muonHitCounts:csc"+intLabel+any));
	  plots["muonHitCountscsc"+intLabel]->Fill(mu.userInt("muonHitCounts:csc"+intLabel));
	  plots["muonHitCountsrpc"+intLabel+any]->Fill(mu.userInt("muonHitCounts:rpc"+intLabel+any));
	  plots["muonHitCountsrpc"+intLabel]->Fill(mu.userInt("muonHitCounts:rpc"+intLabel));
	  rpcTotalHits += mu.userInt("muonHitCounts:rpc"+intLabel);
	}
	
	if(mu.userInt("muonHitCounts") > 0) plots["muonHitCountsratio"]->Fill(float(mu.userInt("muonHitCounts:v1"))/float(mu.userInt("muonHitCounts")));
	if(rpcTotalHits > 0) plots["muonHitCountsrpcratio"]->Fill(float(mu.userInt("muonHitCounts:rpc1"))/float(rpcTotalHits));
	
      }
      
      //for(size_t j =0; j<4; ++j) {
      //std::string intLabel = lexical_cast<std::string>(j+1);
      //  
      //}
      
      float z = 0.;
      float rho = 0.;
      float d = 0.;
      
      if(mu.hasUserFloat("classByHitsGlb:prodZ")) z = mu.userFloat("classByHitsGlb:prodZ");
      if(mu.hasUserFloat("classByHitsGlb:prodRho")) rho = mu.userFloat("classByHitsGlb:prodRho");
      d = sqrt(z*z + rho*rho);
      
      plots["prodz"]->Fill(z);
      plots["prodr"]->Fill(rho);
      plots["prodd"]->Fill(d);
      plots["prodrz"]->Fill(z,rho);
      
    } 
}

void InclusiveMuonPlotsMRTU::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) 
{
  if (luminosity != 0) {
    edm::Handle<edm::MergeableCounter> mc;
    iLumi.getByLabel(normalization_, mc);
    luminosity->Fill(0.5, double(mc->value));
    // set the correct uncertainty from counting statistics
    luminosity->SetBinError(1, sqrt(luminosity->GetBinContent(1)));
  }
}

DEFINE_FWK_MODULE(InclusiveMuonPlotsMRTU);







