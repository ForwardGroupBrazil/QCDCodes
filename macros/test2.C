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
  
  TList * fileList = makeFileCollection("my2112FileList.txt");
  
  fileList->Print();

  TString directories[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/MultiTrack/globalMuons_tpToGlbAssociation"
  }

  TList * dirList = makeDirectoryCollection(fileList,directories,1);

  TIter iter(dirList);
  TDirectory * c;
  while ( (c= (TDirectory *)iter()) ) {
    cout << "The Directory Collection " << c->GetName() << endl;
  }

  //  TList * collection = makeObjectCollection(dirList,"ptres_vs_eta_Sigma");
  TList * collection = makeObjectCollection(dirList,"ptres_vs_eta_Sigma");
  collection->Print();

  drawObjectCollection(collection,false);

  TList * graphList = new TList();

  TIter iterG(collection);
  TH1* h;
  while( (h=(TH1*)iterG()) ) {
    graphList->Add(convertTH1toTGraph(h));
  }
  
  new TCanvas();
  TString legendPt[] = {"muPt10","muPt100","muPt200","muPt500","muPt1000"}
  drawObjectCollection(graphList,true,legendPt);
  ((TGraph*)graphList->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.0001,1.);

  TList * collectionRes = makeObjectCollection(dirList,"ptres_vs_eta");
  TList * resList = makeFitCollection(collectionRes,2);
  ((TGraph*)resList->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.0001,1.);
  TCanvas * canvas = newCanvas("c_test","Test 2 Canvas");

  TLegend * myLegend2 = drawObjectCollection(resList,true,legendPt);

  myLegend2->SetX1(0.2);
  myLegend2->SetX2(0.5);
  myLegend2->SetY1(0.7);
  myLegend2->SetY2(0.9);

  //to delte all items in the list
  //graphList->Delete();

  //to remove but NOT delete all items in the list (unless list owns the objects via list->SetOwner())
  //resList->Clear();

  //to remove the list and all objects in the list
  //resList->SetOwner();
  //delete resList;

  //printCanvasesType(".eps");
  return myLegend2;
  
}
