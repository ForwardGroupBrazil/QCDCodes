TLegend * drawDistributions(TString &dir="./")
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
  TString directories[] = {
    "test_hist/dir1",
    "test_hist/dir2",
    "test_hist/dir6"
  }
  
  TList * fileList = new TList();
  TFile *f;
  f = new TFile(dir,"READ");
  fileList->Add(f);
  fileList->Print();
  TList * dirList = makeDirectoryCollection(fileList,directories,3);
  dirList->Print();

  TString histos[] = {"l3Pt","deltaPt",
		      "glb_calComp","glb_segComp",
		      "glb_nChamber","glb_nChamberMatch",
		      "glb_nCSC_total","glb_nDT_total",
		      "mu_TMLastStationTight","mu_TM2DCompatibilityTight",
		      "mu_GlobalMuonPromptTight",
		      "glb_trkDeltaRelChi2","glb_staDeltaRelChi2",
		      "glb_nSeg_total","glb_nChi2"
  };

  TString legend[] = {"Primary","Silicon","Punch-through"}

  TLegend * theLegend;
  TList * graphCol;
  TList * cutEffCol;

  TCanvas * erm_canvas;
  TCanvas * effGT_canvas;
  TCanvas * effLT_canvas;

  for (int iQuantity=0;iQuantity < 15; ++iQuantity) {
    erm_canvas = newCanvas("/DrafStuffFigs/"+histos[iQuantity],"Name "+histos[iQuantity]);    
    erm_objCol =  makeObjectCollection(dirList,histos[iQuantity]);
    theLegend = drawObjectCollection(erm_objCol,true,legend);
    //graphCol = makeTGfromTHCollection(erm_objCol);
    //theLegend = drawObjectCollection(graphCol,true,legend);
    theLegend->Draw("same");

    effGT_canvas = newCanvas("/DrafStuffFigs/effGT_"+histos[iQuantity],"Name Eff GT "+histos[iQuantity]); 
    cutEffCol = makeCutEffFromTHCollection(erm_objCol,true);
    drawObjectCollection(cutEffCol,true,legend);
    theLegend->Draw("same");

    effLT_canvas = newCanvas("/DrafStuffFigs/effLT_"+histos[iQuantity],"Name Eff LT "+histos[iQuantity]); 
    cutEffCol = makeCutEffFromTHCollection(erm_objCol,false);
    drawObjectCollection(cutEffCol,true,legend);
    theLegend->Draw("same");
  }

  
}
