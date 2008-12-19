#ifndef UserCode_GlobalTruncator_TruncAnalyzer_H
#define UserCode_GlobalTruncator_TruncAnalyzer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "PhysicsTools/RecoAlgos/interface/TrackingParticleSelector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

class DQMStore;
class MonitorElement;
class MuonServiceProxy;
class TrackAssociatorBase;

class TruncAnalyzer : public edm::EDAnalyzer
{
 public:
  TruncAnalyzer(const edm::ParameterSet& pset);
  ~TruncAnalyzer();
  
  virtual void beginJob(const edm::EventSetup& eventSetup);
  virtual void endJob();
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

  //int countStations(const reco::TrackRef&);

 private:


 protected:

  unsigned int verbose_;

  edm::InputTag simLabel_;
  edm::InputTag glbMuLabel_;

  edm::InputTag glbMuAssocLabel_;

  std::vector<edm::InputTag> label;

  std::string outputFileName_;
  std::string subDir_;

  MuonServiceProxy * theMuonService;
  DQMStore * theDQM;

  struct HistoDimensions;
  HistoDimensions * hDim_;

  struct TruncME;
  TruncME * truncMuME_;

  std::vector<TruncME*> V_truncME;

  TrackingParticleSelector tpSelector_;
};

#endif

