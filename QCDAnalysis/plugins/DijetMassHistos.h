#ifndef DijetMassHistos_h
#define DijetMassHistos_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "KKousour/QCDAnalysis/interface/QCDJet.h"
#include "KKousour/QCDAnalysis/interface/QCDEvent.h"
#include "KKousour/QCDAnalysis/interface/QCDEventHdr.h"
#include "KKousour/QCDAnalysis/interface/QCDCaloJet.h"
#include "KKousour/QCDAnalysis/interface/QCDPFJet.h"
#include "KKousour/QCDAnalysis/interface/QCDMET.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"

class DijetMassHistos : public edm::EDAnalyzer 
{
  public:
    explicit DijetMassHistos(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~DijetMassHistos() {}

  private:  
    int getBin(double x, const std::vector<double>& boundaries); 
    //---- configurable parameters --------   
    bool mUsePF;
    double mMinPt1,mMinPt2;
    std::string mFileName,mTreeName;
    std::vector<double> mYBND,mMASSBND,mMinMass,mMaxMass;
    
    edm::Service<TFileService> fs;
    TTree *mTree; 
    TFile *mInf;
    std::vector<TH1F*> mhMETovSUMET,mhM,mhNormM,mhTruncM,mhNormTruncM,mhPt,mhY,mhYmax;
    std::vector<TH1F*> mhCHF,mhNHF,mhPHF,mhN90hits,mhEMF,mhNTrkCalo,mhNTrkVtx,mhfHPD;
    //---- TREE variable --------
    QCDEvent *mEvent;
    
};

#endif
