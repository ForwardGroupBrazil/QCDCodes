#include "UserCode/GlobalTruncator/src/TruncAnalyzer.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"

#include <DataFormats/Common/interface/AssociationMap.h>
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"

#include "PhysicsTools/RecoAlgos/interface/TrackingParticleSelector.h"


#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

using namespace std;
using namespace edm;
using namespace reco;

struct TruncAnalyzer::HistoDimensions {
  unsigned int nBinRes;
  double minResPt, maxResPt;
};

struct TruncAnalyzer::TruncME {
  void bookHistograms(DQMStore* dqm, const string& dirName,  const HistoDimensions& hDim)
  {
    dqm->cd();
    dqm->setCurrentFolder(dirName.c_str());

    hResPtBarrel_  = dqm->book1D("ResPtBarrel" , "#Delta(p_{T})/p_{T}", hDim.nBinRes, hDim.minResPt , hDim.maxResPt );
    hResPtOverlap_  = dqm->book1D("ResPtOverlap" , "#Delta(p_{T})/p_{T}", hDim.nBinRes, hDim.minResPt , hDim.maxResPt );
    hResPtEndcap_  = dqm->book1D("ResPtEndcap" , "#Delta(p_{T})/p_{T}", hDim.nBinRes, hDim.minResPt , hDim.maxResPt );

    hResPBarrel_  = dqm->book1D("ResPBarrel" , "#Delta(p)/p", hDim.nBinRes, hDim.minResPt , hDim.maxResPt );
    hResPOverlap_  = dqm->book1D("ResPOverlap" , "#Delta(p)/p", hDim.nBinRes, hDim.minResPt , hDim.maxResPt );
    hResPEndcap_  = dqm->book1D("ResPEndcap" , "#Delta(p)/p", hDim.nBinRes, hDim.minResPt , hDim.maxResPt );

    hNStations_ = dqm->book1D("NStations","N Stations",10,0,10);
    hLastStation_ = dqm->book1D("LastStation","Last Station",10,0,10);
  };

  //void fill(const TrackingParticle* simRef, const reco::Track* recoRef)
  void fill(const TrackingParticle* simRef, const reco::TrackRef& recoRef)
  {
    const double simP  = simRef->p();
    const double simPt  = simRef->pt();
    const double simEta  = simRef->eta();
    const double recoP  = (recoRef->p());
    const double recoPt  = sqrt(recoRef->momentum().perp2());

    const double errP  = (recoP-simP)/simP;
    const double errPt  = (recoPt-simPt)/simPt;

    if(simEta <= 0.8) hResPtBarrel_ ->Fill(errPt );
    if(simEta > 0.8 && simEta <= 1.2) hResPtOverlap_ ->Fill(errPt );
    if(simEta > 1.2 && simEta <= 2.4) hResPtEndcap_ ->Fill(errPt );

    if(simEta <= 0.8) hResPBarrel_ ->Fill(errP );
    if(simEta > 0.8 && simEta <= 1.2) hResPOverlap_ ->Fill(errP );
    if(simEta > 1.2 && simEta <= 2.4) hResPEndcap_ ->Fill(errP );

    std::pair<int,int> nStation = countStations(recoRef);
    LogDebug("TruncAnalyzer")<<nStation.first<<" " <<nStation.second;
    hNStations_->Fill(nStation.first);
    hLastStation_->Fill(nStation.second);
  };

std::pair<int,int> countStations(const reco::TrackRef& track)
{
  bool dt_1 = false;
  bool dt_2 = false;
  bool dt_3 = false;
  bool dt_4 = false;
  bool csc_1 = false;
  bool csc_2 = false;
  bool csc_3 = false;
  bool csc_4 = false;

  for (trackingRecHit_iterator iall = track->recHitsBegin(); iall != track->recHitsEnd(); ++iall) {
    
    //    if( (*iall)->det()->geographicalId().det()==2 && (*iall)->det()->geographicalId().subdetId()==3) continue;

    if((*iall)->isValid()){
      int DT_station = 0;
      int CSC_station = 0;
      
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
  
  int nTot = (n_dt+n_csc);
 
  return std::pair<int,int>(nTot,last);
};


  typedef MonitorElement* MEP;

  MEP hResPtBarrel_, hResPtOverlap_, hResPtEndcap_, hNStations_, hLastStation_;
  MEP hResPBarrel_, hResPOverlap_, hResPEndcap_;

};

TruncAnalyzer::TruncAnalyzer(const ParameterSet& pset)
{
  verbose_ = pset.getUntrackedParameter<unsigned int>("verbose", 0);

  outputFileName_ = pset.getUntrackedParameter<string>("outputFileName", "");

  minStations_ = pset.getUntrackedParameter<unsigned int>("minStations",4);

  // Set histogram dimensions
  //HistoDimensions hDim;
  hDim_ = new HistoDimensions;
  hDim_->nBinRes  = pset.getUntrackedParameter<unsigned int>("nBinRes");
  hDim_->minResPt = pset.getUntrackedParameter<double>("minResPt");
  hDim_->maxResPt = pset.getUntrackedParameter<double>("maxResPt");

  // Labels for simulation and reconstruction tracks
  simLabel_  = pset.getParameter<InputTag>("simLabel" );
  glbMuLabel_ = pset.getParameter<InputTag>("glbMuLabel");

  // Labels for sim-reco association
  glbMuAssocLabel_ = pset.getParameter<InputTag>("glbMuAssocLabel");

  // the service parameters
  ParameterSet serviceParameters 
    = pset.getParameter<ParameterSet>("ServiceParameters");
  theMuonService = new MuonServiceProxy(serviceParameters);

  label = pset.getParameter< std::vector<edm::InputTag> >("label");


  // retrieve the instance of DQMService
  theDQM = 0;
  theDQM = Service<DQMStore>().operator->();

  if ( ! theDQM ) {
    LogError("TruncAnalyzer") << "DQMService not initialized\n";
    return;
  }

  subDir_ = pset.getUntrackedParameter<string>("subDir");
  if ( subDir_.empty() ) subDir_ = "RecoMuonV";
  //if ( subDir_[subDir_.size()-1] == '/' ) subDir_.erase(subDir_.size()-1);

    tpSelector_ = TrackingParticleSelector(pset.getParameter<double>("ptMinTP"),
					  pset.getParameter<double>("minRapidityTP"),
					  pset.getParameter<double>("maxRapidityTP"),
					  pset.getParameter<double>("tipTP"),
					  pset.getParameter<double>("lipTP"),
					  pset.getParameter<int>("minHitTP"),
					  pset.getParameter<bool>("signalOnlyTP"),
					  pset.getParameter<bool>("chargedOnlyTP"),
					  pset.getParameter<std::vector<int> >("pdgIdTP"));


  if ( verbose_ > 0 ) theDQM->showDirStructure();

}

TruncAnalyzer::~TruncAnalyzer()
{
  if ( theMuonService ) delete theMuonService;
}

void TruncAnalyzer::beginJob(const EventSetup& eventSetup)
{
  if ( theMuonService ) theMuonService->update(eventSetup);
  
  for (unsigned int www=0;www<label.size();www++){
    LogDebug("TruncAnalyzer")<<label[www]; 
      theDQM->cd();
      InputTag algo = label[www];
      string dirName=subDir_;
      if (algo.process()!="")
	dirName+=algo.process()+"_";
      if(algo.label()!="")
	dirName+=algo.label()+"_";
      if(algo.instance()!="")
	dirName+=algo.instance()+"_";      

      std::replace(dirName.begin(), dirName.end(), ':', '_');
      //      theDQM->setCurrentFolder(dirName.c_str());

      theDQM->cd();
      theDQM->setCurrentFolder(dirName.c_str());

      TruncME * tme = new TruncME;
      tme->bookHistograms(theDQM,dirName.c_str(),*hDim_);
      V_truncME.push_back(tme);

 }

  //glbMuAssociator_ = 0;

}

void TruncAnalyzer::endJob()
{
  if ( theDQM && ! outputFileName_.empty() ) theDQM->save(outputFileName_);
}

void TruncAnalyzer::analyze(const Event& event, const EventSetup& eventSetup)
{
  using namespace reco;
  // Get TrackingParticles
  Handle<TrackingParticleCollection> simHandle;
  event.getByLabel(simLabel_, simHandle);
  const TrackingParticleCollection simColl = *(simHandle.product());

  // Get Muon Tracks
  Handle<View<Track> > glbColl;
  event.getByLabel(glbMuLabel_, glbColl);
  //const View<Track> glbMuColl = *(glbMuHandle.product());
  
  // Get Association maps
  SimToRecoCollection simToGlbColl;
  RecoToSimCollection glbToSimColl;
  
  Handle<SimToRecoCollection> simToGlbMuHandle;
  event.getByLabel(glbMuAssocLabel_, simToGlbMuHandle);
  simToGlbColl = *(simToGlbMuHandle.product());
  
  Handle<RecoToSimCollection> glbMuToSimHandle;
  event.getByLabel(glbMuAssocLabel_, glbMuToSimHandle);
  glbToSimColl = *(glbMuToSimHandle.product());
  LogDebug("TruncAnalyzer") << simColl.size();
  LogDebug("TruncAnalyzer") << simToGlbColl.size();
  for(TrackingParticleCollection::size_type i=0; i<simColl.size(); i++) {
    TrackingParticleRef simRef(simHandle, i);
    const TrackingParticle* simTP = simRef.get();

    if ( ! tpSelector_(*simTP) ) continue;
    //const double simPt  = simRef->pt();

    // Get sim-reco association for a simRef
    vector<pair<RefToBase<Track>, double> > glbMuRefV;

    if ( simToGlbColl.find(simRef) != simToGlbColl.end() ) {

      glbMuRefV = simToGlbColl[simRef];

      if ( ! glbMuRefV.empty() ) {
	const reco::TrackRef glbMuTrackRef = glbMuRefV.begin()->first.castTo<TrackRef >();

	//
	std::pair<int,int> count = V_truncME[0]->countStations(glbMuTrackRef);
	//

	if(count.first >= minStations_)	
	  for (unsigned int www=0;www<label.size();www++){
	    
	    edm::Handle<reco::TrackToTrackMap> truncAssoMap;
	    bool map = event.getByLabel(label[www],truncAssoMap);
	    //LogDebug("TruncAnalyzer")<<label[www]<<" " <<map;
	    reco::TrackToTrackMap::const_iterator iEnd;
	    reco::TrackToTrackMap::const_iterator iii;
	    if(map) {
	      iEnd = truncAssoMap->end();
	      iii = truncAssoMap->find(glbMuTrackRef);
	      //const Track * tk = (*truncAssoMap)[glbMuTrackRef].get();
	      const TrackRef  tk = (*truncAssoMap)[glbMuTrackRef];
	      V_truncME[www]->fill(simTP,tk);
	    }	  	  	  
	  }
      }
    }
  }
}
