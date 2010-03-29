TLegend * drawTFS(TString &dir="./")
{
  //----------------------------------------------------
  // General macro setup
  //----------------------------------------------------
  
  gROOT->LoadMacro("adamStyles.C");
  
  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();
  
  gROOT->LoadMacro("adamCanvasMacros.C");
  gROOT->LoadMacro("adamDrawMacros.C");
  gROOT->LoadMacro("adamFitMacros.C");
  gROOT->LoadMacro("adamGetObjMacros.C");
  gROOT->LoadMacro("adamMakeCollectionMacros.C");

  //----------------------------------------------------
  // General setup for files and directories
  //----------------------------------------------------

  TList * fileList = makeFileCollection("my09BadFileList.txt");
  
  fileList->Print();

  TString directories[] = {
    //"globalMuons/trackMatcher"
    //"globalMuons/builderBase"
    //"globalMuons/builder"
    "standAloneMuons/trackLoader"
    //"globalMuons/trackLoader"
  }
  
  //TList * fileList = new TList();
  //TFile *f;
  //f = new TFile(dir,"READ");
  //fileList->Add(f);
  //fileList->Print();
  TList * dirList = makeDirectoryCollection(fileList,directories,1);
  dirList->Print();

//   TString histos[] = {"h_pt",
// 		      "h_nTk",
// 		      "h_distance_0",
// 		      "h_chi2_0",
// 		      "h_loc_chi2_0",
// 		      "h_deltaR_0",
// 		      "h_distance_a1",
// 		      "h_loc_chi2_a1",
// 		      "h_chi2_a2",
// 		      "h_distance_a2",
// 		      "h_distance_a3",
// 		      "h_deltaR_a3",
// 		      "h_deltaEta_b",
// 		      "h_deltaPhi_b",
// 		      "h_nMatch",
// 		      "h_nClean",
// 		      "h_nCands",
// 		      "h_na1",
// 		      "h_na2",
// 		      "h_na3",
// 		      "h_nb"
//   } //21

//   TString histos[] = {
//     "h_nTkTrajs",
//     "h_nStaTkRefittedTrajs",
//     "h_staTkProb",
//   } // 3

//   TString histos[] = {
//     "h_nRegionalTk",
//     "h_nMatchedTk",
//     "h_nSta",
//     "h_nGlb",
//     "h_staPt",
//     "h_staRho",
//     "h_staR",
//   } //7

  TString histos[] = {
    "h_dInner",
    //    "h_dPerp"
  } //1


  TString legend[] = {"Data","MC"}

  TLegend * theLegend;
  TList * graphCol;
  TList * cutEffCol;

  TCanvas * erm_canvas;
  TCanvas * effGT_canvas;
  TCanvas * effLT_canvas;

  for (int iQuantity=0;iQuantity < 1; ++iQuantity) {
    erm_canvas = newCanvas("./DrawStuffFigs_Bad/"+directories[0]+"/"+histos[iQuantity],"Name "+histos[iQuantity]);    
    erm_objCol =  makeObjectCollection(dirList,histos[iQuantity]);
    theLegend = drawObjectCollection(erm_objCol,true,legend);
    theLegend->SetX1(0.7);
    theLegend->SetX2(0.9);
    theLegend->SetY1(0.65);
    theLegend->SetY2(0.85);
    //graphCol = makeTGfromTHCollection(erm_objCol);
    //theLegend = drawObjectCollection(graphCol,true,legend);
    erm_canvas->SetLogy(1);
    theLegend->Draw("same");

//     effGT_canvas = newCanvas("./DrawStuffFigs/effGT_"+histos[iQuantity],"Name Eff GT "+histos[iQuantity]); 
//     cutEffCol = makeCutEffFromTHCollection(erm_objCol,true);
//     drawObjectCollection(cutEffCol,true,legend);
//     theLegend->Draw("same");

//     effLT_canvas = newCanvas("./DrawStuffFigs/effLT_"+histos[iQuantity],"Name Eff LT "+histos[iQuantity]); 
//     cutEffCol = makeCutEffFromTHCollection(erm_objCol,false);
//     drawObjectCollection(cutEffCol,true,legend);
//     theLegend->Draw("same");
  }

  //  printCanvasesType(".eps");
}
