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
    bool mIsMC;
    double mMinPt;
    std::string mFileName,mTreeName,mDirName;
    std::vector<double> mYBND,mPTBND;
    std::vector<std::string> mTriggers;
    std::vector<int> mTrigIndex;    

    edm::Service<TFileService> fs;
    TTree *mTree; 
    TFile *mInf;
    TDirectoryFile *mDir;
    TH1F *mhJetMulti[30];
    TH1F *mhPt[30][10],*mhNormPt[30][10],*mhX[30][10],*mhNormX[30][10],
         *mhCHF[30][10],*mhNHF[30][10],*mhPHF[30][10],
         *mhN90hits[30][10],*mhEMF[30][10],*mhNTrkCalo[30][10],*mhNTrkVtx[30][10],*mhfHPD[30][10],*mhMETovSUMET[30][10];
    //---- TREE variable --------
    QCDEvent *mEvent;
    
};

#endif
