#include "FillHistograms.C"
void FillHistogramsAll()
{
  TString FILENAME[12] = {
    "/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_data_preselect_hard_tmva",
    "/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_QCD-HT100_preselect_hard_tmva",
    "/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_QCD-HT250_preselect_hard_tmva",
    "/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_QCD-HT500_preselect_hard_tmva",
    "/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_QCD-HT1000_preselect_hard_tmva",
    "/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_ZJets_preselect_hard_tmva",
    "/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_TTJets_preselect_hard_tmva",
    "/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_T_preselect_hard_tmva",
    "/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_Tbar_preselect_hard_tmva",
    "/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_VBF-Powheg125_preselect_hard_tmva",
    "/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_GluGlu-Powheg125_preselect_hard_tmva",
    "/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_GluGlu-Madgraph125_preselect_hard_tmva"
  };
  bool ApplyTriggerSel[12] = {true,true,true,true,true,true,true,true,true,true,true,true};
  bool isMC[12]            = {false,true,true,true,true,true,true,true,true,true,true,true}; 
  
  for(int i=0;i<12;i++) {
    FillHistograms(FILENAME[i],ApplyTriggerSel[i],isMC[i]);
  }  
}
