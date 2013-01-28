#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TKey.h"
#include "TCollection.h"
#include "TDirectoryFile.h"
#include <iostream>
using std::cin;
using std::cout;
using std::endl;
void PreselectSamples(TString PRESELECTION)
{
  TString PRE("");
  if (PRESELECTION.CompareTo("soft") == 0) {
    PRE = "jetPt[0]>85 && jetPt[1]>70 && jetPt[2]>50 && jetPt[3]>35 && dEtaqq>2.5 && mqq>250";
  }
  else if (PRESELECTION.CompareTo("hard") == 0) {
    PRE = "jetPt[0]>85 && jetPt[1]>70 && jetPt[2]>60 && jetPt[3]>40 && dEtaqq>2.5 && mqq>300";
  }
  else {
    cout<<"Wrong preselection type !!! Must be \"soft\" or \"hard\" "<<endl;
    return;
  }
  TString FileName[19] = {
  "flatTree_MultiJet-Run2012A","flatTree_BJetPlusX-Run2012B","flatTree_BJetPlusX-Run2012C","flatTree_BJetPlusX-Run2012D",
  "flatTree_VBF-Powheg115","flatTree_VBF-Powheg120","flatTree_VBF-Powheg125","flatTree_VBF-Powheg130","flatTree_VBF-Powheg135",
  "flatTree_GluGlu-Powheg125","flatTree_GluGlu-Madgraph125",
  "flatTree_QCD-HT100","flatTree_QCD-HT250","flatTree_QCD-HT500","flatTree_QCD-HT1000",
  "flatTree_Zjets","flatTree_T","flatTree_Tbar","flatTree_TTJets"
  };
  TString VAR[19] = {"pvx","pvy","pvz","ht","jetRChg_QC","softTrackJetPt","softTrackJetEta","softTrackJetPhi","softTrackJetE",
  "genjetPt","genjetEta","genjetPhi","genjetE","partonId","partonSt","partonPt","partonEta","partonPhi","partonE"};
  TKey *key;
  for(int iFile=0;iFile<15;iFile++) {
    cout<<"Opening file: "<<FileName[iFile]<<endl;
    TFile *inf = TFile::Open("root://eoscms//eos/cms/store/cmst3/user/kkousour/VBF-CMG/"+FileName[iFile]+".root");
    // find the directories
    std::vector<TString> DIR_NAME;
    TIter nextkey(inf->GetListOfKeys());
    while ((key = (TKey*)nextkey())) {
      if (TString(key->GetClassName()).CompareTo("TDirectoryFile") == 0) {
        DIR_NAME.push_back(TString(key->GetName()));
      }  
    }
    TFile *outf  = TFile::Open(FileName[iFile]+"_preselect_"+PRESELECTION+".root","RECREATE");
    for(unsigned int idir=0;idir<DIR_NAME.size();idir++) {
      TTree *tIN = (TTree*)inf->Get(DIR_NAME[idir]+"/events");
      //----- disable variables ----------------
      for(int ivar=0;ivar<19;ivar++) {
        tIN->SetBranchStatus(VAR[ivar],0);
      }      
      cout<<"Cloning tree: "<<DIR_NAME[idir]<<" with "<<tIN->GetEntries()<<" entries"<<endl;
      TTree *trOUT = (TTree*)tIN->CopyTree(PRE.Data());
      delete tIN;
      TDirectoryFile *dir = (TDirectoryFile*)outf->mkdir(DIR_NAME[idir]);
      dir->cd();
      trOUT->Write("events");
      cout<<"Wrote out tree with "<<trOUT->GetEntries()<<" entries"<<endl;
      dir->Close();
      delete dir;
      delete trOUT;
    }  
    outf->Close();
    inf->Close();
  }// end of file loop
}
