#ifndef RecoMuon_StandAloneMuonProducer_STAMuonAnalyzer2_H
#define RecoMuon_StandAloneMuonProducer_STAMuonAnalyzer2_H

/** \class STAMuonAnalyzer2
 *  Analyzer of the StandAlone muon tracks
 *
 *  $Date: 2006/09/01 14:35:48 $
 *  $Revision: 1.2 $
 *  \author R. Bellan - INFN Torino <riccardo.bellan@cern.ch>
 */

// Base Class Headers
#include "FWCore/Framework/interface/EDAnalyzer.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

class TFile;
class TH1F;
class TH2F;

class STAMuonAnalyzer2: public edm::EDAnalyzer {
public:
  /// Constructor
  STAMuonAnalyzer2(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~STAMuonAnalyzer2();

  // Operations

  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);

  virtual void beginJob(const edm::EventSetup& eventSetup) ;
  virtual void endJob() ;
protected:

private:
  std::string theRootFileName;
  TFile* theFile;

  std::string theSTAMuonLabel;
  std::string theSeedCollectionLabel;

  // Histograms
  TH1F *hPtRec;
  TH1F *hPtSim; 
  TH1F *hPres;
  TH1F *h1_Pres;
  TH1F *hPTDiff;
  TH1F *hPTDiff2;
  TH2F *hPTDiffvsEta;
  TH2F *hPTDiffvsPhi;
  
  TH1F *hNtracks;
  TH1F *hNsim;
  TH1F *hNrec;
  TH1F *hNmatch;
  TH1F *hetasim;
  TH1F *hptsim;
  TH1F *hdR;
  TH1F *hetarec;
  TH1F *hetarec1;
  TH1F *hetarec_all;
  TH1F *hptrec;
  TH1F *hptrec1;
  TH2F *hdptrecvsEta;
  TH2F *hptrecvsEta;
  TH1F *hdptrec;
  TH1F *hminv; 
  TH1F *hminvrec; 
  TH1F *hminvrecZ; 

  TH1F* hetaeff;
  TH1F* hetapur;

  // Counters
  int numberOfSimTracks;
  int numberOfRecTracks;

  std::string theDataType;
  
};
#endif

