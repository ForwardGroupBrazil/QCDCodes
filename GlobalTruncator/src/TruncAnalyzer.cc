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

    hResPt_  = dqm->book1D("ResPt" , "#Delta(p_{T})/p_{T}", hDim.nBinRes, hDim.minResPt , hDim.maxResPt );
    hNStations_ = dqm->book1D("NStations","N Stations",10,0,10);
  };

  void fill(const TrackingParticle* simRef, const Track* recoRef)
  {
    const double simPt  = simRef->pt();
    const double recoPt  = sqrt(recoRef->momentum().perp2());

    const double errPt  = (recoPt-simPt)/simPt;

    hResPt_ ->Fill(errPt );

    //hNStations_->Fill(countStations(recoRef));
  };

  typedef MonitorElement* MEP;

  MEP hResPt_, hNStations_;

};

TruncAnalyzer::TruncAnalyzer(const ParameterSet& pset)
{
  verbose_ = pset.getUntrackedParameter<unsigned int>("verbose", 0);

  outputFileName_ = pset.getUntrackedParameter<string>("outputFileName", "");

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
	TrackRef glbMuTrackRef = glbMuRefV.begin()->first.castTo<TrackRef >();
	
	for (unsigned int www=0;www<label.size();www++){
	  
	  edm::Handle<reco::TrackToTrackMap> truncAssoMap;
	  bool map = event.getByLabel(label[www],truncAssoMap);

	  reco::TrackToTrackMap::const_iterator iEnd;
	  reco::TrackToTrackMap::const_iterator iii;
	  if(map) {
	    iEnd = truncAssoMap->end();
	    iii = truncAssoMap->find(glbMuTrackRef);
	    const Track * tk = (*truncAssoMap)[glbMuTrackRef].get();
	    V_truncME[www]->fill(simTP,tk);
	  }
	  

	  
	}
      }
    }
    
    }
    
    
}


int
TruncAnalyzer::countStations(const reco::Track* track)
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
  int n_dt = 0;
  if (dt_1) n_dt++;
  if (dt_2) n_dt++;
  if (dt_3) n_dt++;
  if (dt_4) n_dt++;

  int n_csc = 0;
  if (csc_1) n_csc++;
  if (csc_2) n_csc++;
  if (csc_3) n_csc++;
  if (csc_4) n_csc++;
  
  return ((n_dt+n_csc) == 1 && (dt_1 || csc_1));
}




