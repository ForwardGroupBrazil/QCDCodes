#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TMVA/Factory.h"
using namespace TMVA;
void RegressionMVA()
{
  TFile *inf = TFile::Open("SignalRegressionTrain.root");
  TTree *tr  = (TTree*)inf->Get("jets");
  
  TFile *outf = new TFile("regressionTMVA.root","RECREATE");
  
  TMVA::Factory* factory = new TMVA::Factory("factoryJetReg",outf,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression");

  factory->AddRegressionTree(tr);
  
  factory->AddVariable("jetBtag"   ,'F');
  //factory->AddVariable("jetE"      ,'F');
  factory->AddVariable("jetPt"     ,'F');
  factory->AddVariable("jetEta"    ,'F');
  factory->AddVariable("jetMetPhi" ,'F');
  factory->AddVariable("jetChf"    ,'F');
  factory->AddVariable("jetPhf"    ,'F');
  factory->AddVariable("jetNhf"    ,'F');
  factory->AddVariable("jetElf"    ,'F');
  factory->AddVariable("jetMuf"    ,'F');
  factory->AddVariable("jetPtD"    ,'F');
  factory->AddVariable("jetVtxPt"  ,'F');
  factory->AddVariable("jetVtx3dL" ,'F');
  factory->AddVariable("jetVtx3deL",'F');
  factory->AddVariable("met"       ,'F');
  factory->AddVariable("rho"       ,'F');
  
  factory->AddTarget("jetGenPt");
    
  TCut preselectionCut("jetGenDR<0.25");

  factory->PrepareTrainingAndTestTree(preselectionCut,"nTrain_Regression=60000:nTest_Regression=60000");
  factory->BookMethod(TMVA::Types::kMLP,"MLP","NCycles=700:HiddenLayers=N,N-1:TestRate=5:TrainingMethod=BFGS:VarTRansform=Norm");
  factory->BookMethod(TMVA::Types::kBDT,"BDT","NTrees=200;nCuts=25"); 
    
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods(); 
}
