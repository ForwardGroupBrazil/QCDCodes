TLegend * test2()
{
  gROOT->LoadMacro("adamStyles.C");
  
  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();
  
  gROOT->LoadMacro("adamCanvasMacros.C");
  gROOT->LoadMacro("adamDrawMacros.C");
  gROOT->LoadMacro("adamFitMacros.C");
  gROOT->LoadMacro("adamGetObjMacros.C");
  gROOT->LoadMacro("adamMakeCollectionMacros.C");
  
  TSortedList * fileList = makeFileCollection("myTestFileList.txt");
  
  fileList->Print();

  TString directories[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/MultiTrack/general_tpToTkmuAssociation",
    "/DQMData/Run 1/RecoMuonV/Run summary/MultiTrack/globalMuons_tpToGlbAssociation"
  }

  TSortedList * dirList = makeDirectoryCollection(fileList,directories,2);

  TIter iter(dirList);
  TDirectory * c;
  while ( (c= (TDirectory *)iter()) ) {
    cout << "The Directory Collection " << c->GetName() << endl;
  }

  //  TSortedList * collection = makeObjectCollection(dirList,"ptres_vs_eta_Sigma");
  TSortedList * collection = makeObjectCollection(dirList,"effic");
  collection->Print();

  drawObjectCollection(collection,false);

  TSortedList * graphList = new TSortedList();

  TIter iterG(collection);
  TH1* h1;
  while( (h=(TH1*)iterG()) ) {
    graphList->Add(convertTH1toTGraph(h));
  }

  new TCanvas();
  drawObjectCollection(graphList);

  TSortedList * collectionRes = makeObjectCollection(dirList,"ptres_vs_eta");
  TSortedList * resList = makeFitCollection(collectionRes,2);
  ((TGraph*)resList->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.01,0.08);
  TCanvas * canvas = newCanvas("c_test","Test 2 Canvas");
  TString legend[] = {"General Tracks","Global Muons"}
  TLegend * myLegend2 = drawObjectCollection(resList,true,legend);

  //printCanvasesType(".eps");
  return myLegend2;
}
