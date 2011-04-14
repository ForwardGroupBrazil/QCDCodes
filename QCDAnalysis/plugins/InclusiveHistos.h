#ifndef InclusiveHistos_h
#define InclusiveHistos_h

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

class InclusiveHistos : public edm::EDAnalyzer 
{
  public:
    explicit InclusiveHistos(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~InclusiveHistos();

  private:  
    int getBin(double x, const std::vector<double>& boundaries); 
    //---- configurable parameters --------   
    bool mUsePF;
    std::string mFileName,mTreeName;
    std::vector<double> mYBND,mPTBND,mMinPt,mMaxPt;
    
    edm::Service<TFileService> fs;
    TTree *mTree; 
    TFile *mInf;
    TH1F *mhMETovSUMET;
    std::vector<TH1F*> mhPt,mhNormPt,mhTruncPt,mhNormTruncPt,mhCHF,mhNHF,mhPHF,mhN90hits,mhEMF,mhNTrkCalo,mhNTrkVtx,mhfHPD;
    //---- TREE variable --------
    QCDEvent *mEvent;
    
};

#endif
