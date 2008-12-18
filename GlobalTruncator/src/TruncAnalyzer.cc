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

using namespace std;
using namespace edm;
using namespace reco;

struct HistoDimensions {
  unsigned int nBinRes;
  double minResPt, maxResPt;
};

struct TruncAnalyzer::TruncME {
  void bookHistograms(DQMStore* dqm, const string& dirName,  const HistoDimensions& hDim)
  {
    dqm->cd();
    dqm->setCurrentFolder(dirName.c_str());

    hResPt_  = dqm->book1D("ResPt" , "#Delta(p_{T})/p_{T}", hDim.nBinRes, hDim.minResPt , hDim.maxResPt );
  };

  void fill(const TrackingParticle* simRef, const Track* recoRef)
  {
    const double simPt  = simRef->pt();
    const double recoPt  = sqrt(recoRef->momentum().perp2());

    const double errPt  = (recoPt-simPt)/simPt;

    hResPt_ ->Fill(errPt );
  };

  typedef MonitorElement* MEP;

  MEP hResPt_;

};

TruncAnalyzer::TruncAnalyzer(const ParameterSet& pset)
{
  verbose_ = pset.getUntrackedParameter<unsigned int>("verbose", 0);

  outputFileName_ = pset.getUntrackedParameter<string>("outputFileName", "");

  // Set histogram dimensions
  HistoDimensions hDim;

  hDim.nBinRes  = pset.getUntrackedParameter<unsigned int>("nBinRes");
  hDim.minResPt = pset.getUntrackedParameter<double>("minResPt");
  hDim.maxResPt = pset.getUntrackedParameter<double>("maxResPt");

  // Labels for simulation and reconstruction tracks
  simLabel_  = pset.getParameter<InputTag>("simLabel" );
  glbMuLabel_ = pset.getParameter<InputTag>("glbMuLabel");

  // Labels for sim-reco association
  glbMuAssocLabel_ = pset.getParameter<InputTag>("glbMuAssocLabel");

  // the service parameters
  ParameterSet serviceParameters 
    = pset.getParameter<ParameterSet>("ServiceParameters");
  theMuonService = new MuonServiceProxy(serviceParameters);

  // retrieve the instance of DQMService
  theDQM = 0;
  theDQM = Service<DQMStore>().operator->();

  if ( ! theDQM ) {
    LogError("TruncAnalyzer") << "DQMService not initialized\n";
    return;
  }

  subDir_ = pset.getUntrackedParameter<string>("subDir");
  if ( subDir_.empty() ) subDir_ = "RecoMuonV";
  if ( subDir_[subDir_.size()-1] == '/' ) subDir_.erase(subDir_.size()-1);

  // book histograms
  theDQM->cd();

  theDQM->setCurrentFolder(subDir_+"/Muons");

  //begin adam questions
  truncMuME_ = new TruncME;

  theDQM->setCurrentFolder(subDir_+"/Glb");
  theDQM->bookString("TrackLabel", glbMuLabel_.label()+"_"+glbMuLabel_.instance());
  theDQM->bookString("AssocLabel", glbMuAssocLabel_.label());

  truncMuME_->bookHistograms(theDQM, subDir_+"/Glb", hDim);



  //end adam questions

  if ( verbose_ > 0 ) theDQM->showDirStructure();

}

TruncAnalyzer::~TruncAnalyzer()
{
  if ( theMuonService ) delete theMuonService;
}

void TruncAnalyzer::beginJob(const EventSetup& eventSetup)
{
  if ( theMuonService ) theMuonService->update(eventSetup);

  //glbMuAssociator_ = 0;

}

void TruncAnalyzer::endJob()
{
  if ( theDQM && ! outputFileName_.empty() ) theDQM->save(outputFileName_);
}

void TruncAnalyzer::analyze(const Event& event, const EventSetup& eventSetup)
{

}






