void trainTMVA(TString SIGNAL)
{
  TFile *bkgSrc  = TFile::Open("data/flatTree_qcd_weights.root");
  TFile *sigSrc  = TFile::Open("data/flatTree_"+SIGNAL+"_weights.root");
  TTree* bkgTree = (TTree*)bkgSrc->Get("multijets/events");
  TTree* sigTree = (TTree*)sigSrc->Get("multijets/events"); 
  TFile *outf    = new TFile("data/TMVA_"+SIGNAL+".root","RECREATE");

  TCut preselectionCut = "ht>850 && pt[7]>30";
  int Nbkg = bkgTree->GetEntries();
  int Nsig = sigTree->GetEntries();

  TMVA::Factory* factory = new TMVA::Factory("factory_"+SIGNAL,outf,"!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification" );

  factory->AddSignalTree(sigTree,1.0);
  factory->AddBackgroundTree(bkgTree,1.0);

  factory->SetWeightExpression("wt");
  
  factory->AddVariable("pt0 := pt[0]",'F');
  factory->AddVariable("pt1 := pt[1]",'F');
  factory->AddVariable("pt2 := pt[2]",'F');
  factory->AddVariable("pt3 := pt[3]",'F');
  factory->AddVariable("pt4 := pt[4]",'F');
  factory->AddVariable("pt5 := pt[5]",'F');
  factory->AddVariable("pt6 := pt[6]",'F');
  factory->AddVariable("pt7 := pt[7]",'F');
  factory->AddVariable("b4j := m4jBalance",'F');
  factory->AddVariable("b2j := m2jSig/m2jAve",'F');
  factory->AddVariable("m2j := m2jAve",'F');  
  factory->AddVariable("m4j := m4jAve",'F'); 
  factory->AddVariable("m8j := m8j",'F');
  factory->AddVariable("ht  := ht",'F');
  
  char name[1000];
  sprintf(name,"nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d",Nsig/2.,Nbkg/2.,Nsig/2.,Nbkg/2.);
  
  factory->PrepareTrainingAndTestTree(preselectionCut,name);
  
  //factory->BookMethod(TMVA::Types::kLikelihood,"LikelihoodD","!H:!V:!TransformOutput:PDFInterpol=Spline2:\
NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:\
NAvEvtPerBin=50:VarTransform=Decorrelate:CreateMVAPdfs=True");
  factory->BookMethod(TMVA::Types::kFisher,"Fisher");
  factory->BookMethod(TMVA::Types::kBDT,"BDT");
  factory->BookMethod(TMVA::Types::kMLP,"MLP_ANN","NCycles=200:HiddenLayers=N,N-1:TestRate=5:TrainingMethod=BFGS:VarTRansform=Norm");

  factory->TrainAllMethods();

  factory->TestAllMethods();
  
  factory->EvaluateAllMethods(); 

}
