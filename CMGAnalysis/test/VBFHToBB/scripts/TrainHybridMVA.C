#include "TMVA/Factory.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h" 
using namespace TMVA;
using namespace TMath;
void TrainHybridMVA()
{
  // the training is done using a dedicated tree format
  TFile *bkgSrc  = TFile::Open("/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_data_PF_CALO_reduced_train.root");
  if (bkgSrc->IsZombie()) {
    cout<<"Bkg source file does not exist !!!"<<endl;
    return;
  }
  TFile *sigSrc  = TFile::Open("/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/flatTree_Powheg125_train.root");
  if (sigSrc->IsZombie()) {
    cout<<"Signal source file does not exist !!!"<<endl;
    return;
  }
  TTree* bkgTree = (TTree*)bkgSrc->Get("events");
  TTree* sigTree = (TTree*)sigSrc->Get("events"); 
  TFile *outf    = new TFile("hybridTMVA20_125.root","RECREATE");

  TCut preselectionCut = "jetPt[0]>85 && jetPt[1]>70 && jetPt[2]>60 && jetPt[3]>40 && dEtaqq>2.5 && mbb>50 && mbb<200 && mqq>250 && jetBtag[0]>0 && jetBtag[1]>0 && dPhibb<2.0";
  
  TMVA::Factory* factory = new TMVA::Factory("factory_hybrid20_125",outf,"!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification" );
  
  factory->AddSignalTree(sigTree);
  factory->AddBackgroundTree(bkgTree);
  
  factory->AddVariable("dEtaqq"           ,'F');
  factory->AddVariable("dEtaqqEta-dEtaqq" ,'F');
  factory->AddVariable("mqq"              ,'F');
  factory->AddVariable("etaBoostqq"       ,'F');
  factory->AddVariable("jetBtag[0]"       ,'F');
  factory->AddVariable("jetBtag[1]"       ,'F');
  factory->AddVariable("jetQGLnew[2]"     ,'F');
  factory->AddVariable("jetQGLnew[3]"     ,'F');
  factory->AddVariable("softHt"           ,'F');
  factory->AddVariable("softMulti"        ,'I');
  factory->AddVariable("jetEtaBtag[2]"    ,'F');
  factory->AddVariable("jetPt[4]/mqq"     ,'F');
  // spectator variables: not used for the training but recorded
  factory->AddSpectator("mbb"             ,'F');
  factory->AddSpectator("dPhibb"          ,'F');
  factory->AddSpectator("rho"             ,'F');

  // the number of background events is chosen such that the
  // events left for training are roughly equal to the signal ones
  factory->PrepareTrainingAndTestTree(preselectionCut,"nTrain_Signal=0:nTrain_Background=136000:nTest_Signal=0:nTest_Background=136000");

  // specify the training methods
  factory->BookMethod(TMVA::Types::kFisher,"Fisher");
  //factory->BookMethod(TMVA::Types::kMLP,"MLP_ANN","NCycles=1400:HiddenLayers=N,N-1:TestRate=5:TrainingMethod=BFGS:VarTRansform=Norm");
  factory->BookMethod(TMVA::Types::kBDT,"BDT","NTrees=400:nCuts=25");
  
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods(); 
}
