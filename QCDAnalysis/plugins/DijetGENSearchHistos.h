#ifndef DijetGENSearchHistos_h
#define DijetGENSearchHistos_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"

class DijetGENSearchHistos : public edm::EDAnalyzer 
{
  public:
    explicit DijetGENSearchHistos(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~DijetGENSearchHistos();

  private:  
    //---- configurable parameters --------   
    double mMinPt;
    double mMinMass;
    double mMaxEta;
    double mMaxDeta;
    double mChiIN;
    double mChiOUT; 
    std::string mGenJetsName;
    std::vector<double> mMassBND;

    edm::Service<TFileService> fs;
    TH1F *mhPt,*mhPtJJ,*mhMass,*mhEta,*mhEtaBoost,*mhDeta,*mhDphi,*mhMassIN,*mhMassOUT,*mhChi,*mhChiIN,*mhChiOUT;
    
};

#endif
