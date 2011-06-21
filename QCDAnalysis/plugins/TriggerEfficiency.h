#ifndef TriggerEfficiency_h
#define TriggerEfficiency_h

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

class TriggerEfficiency : public edm::EDAnalyzer 
{
  public:
    explicit TriggerEfficiency(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~TriggerEfficiency();

  private:  
    int getBin(double x, const std::vector<double>& boundaries); 
    //---- configurable parameters --------   
    int mJetID,mHCALNoise,mNEvents;   
    double mMinPt;
    std::vector<int> mRefTrigIndex;
    std::vector<double> mYBND,mL1Pt,mHLTPt;
    std::vector<std::string> mRefTrigger; 
    std::string mFileName,mTreeName,mDirName;

    edm::Service<TFileService> fs;
    TTree *mTree; 
    TFile *mInf;
    TDirectoryFile *mDir;
    TH1F *mPFPt[30][10],*mPFRefPt[30][10],*mCaloPt[30][10],*mCaloRefPt[30][10];
    //---- TREE variable --------
    QCDEvent *mEvent;
    
};

#endif
