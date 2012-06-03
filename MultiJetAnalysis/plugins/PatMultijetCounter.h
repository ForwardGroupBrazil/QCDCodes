#ifndef PatMultijetCounter_h
#define PatMultijetCounter_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TH1F.h"

class PatMultijetCounter : public edm::EDAnalyzer 
{
  public:
    explicit PatMultijetCounter(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~PatMultijetCounter();

  private:  
    //---- configurable parameters --------   
    bool isMC_;
    double etaMax_,ptMin_,htMin_,betaMin_;
    edm::InputTag srcJets_;
    edm::Service<TFileService> fs_;
    TH1F *hJetMulti_; 
};

#endif
