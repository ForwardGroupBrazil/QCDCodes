#ifndef ResponseHistos_h
#define ResponseHistos_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "KKousour/QCDAnalysis/interface/QCDEvent.h"
#include "TTree.h"
#include "TH2F.h"
#include "TFile.h"

class ResponseHistos : public edm::EDAnalyzer 
{
  public:
    explicit ResponseHistos(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~ResponseHistos();

  private:  
    int getBin(double x, const std::vector<double>& boundaries);
    //---- configurable parameters --------   
    std::string mFileName,mTreeName,mDirName;
    std::vector<double> mPtBND;
    std::vector<double> mEtaBND;
    double mMaxDR;
    int    mNEvents;

    edm::Service<TFileService> fs;
    QCDEvent *mEvent;
    TTree *mTree; 
    TFile *mInf;
    TDirectoryFile *mDir;
    TH2F *mCaloRspVsEta[20],*mCaloRspVsPt[20],*mPFRspVsEta[20],*mPFRspVsPt[20];
};

#endif
