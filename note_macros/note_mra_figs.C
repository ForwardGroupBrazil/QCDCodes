TLegend * note_mra_figs()
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

  TString dirName_("");

  TList * fileList = makeFileCollection("MyBinnedList.txt");

  TString directories[] = {
    "/DQMData/Run 1/Muons/Run summary/MuonRecoAnalyzer/"+dirName_
  }

  TList * dirList = makeDirectoryCollection(fileList,directories,1);


  TString legendPt[] = {"muPt100","muPt1000"}

  //----------------------------------------------------
  // Now do each specific figure
  //----------------------------------------------------

  TString level = "GlbMuon_";
  TString missingHistoNames[] = {
    //"Pull_TkSta_eta",
    //"Pull_TkSta_phi",
    //"Pull_TkSta_qOverp",
    //"Pull_TkSta_oneOverp",
    //"Pull_TkSta_qOverpt",
    //"Pull_TkSta_oneOverpt",
  }
  TString histoName[] = {
    "muReco",
    "Res_TkGlb_eta","Res_GlbSta_eta","Res_TkSta_eta",
    "Res_TkGlb_phi","Res_GlbSta_phi","Res_TkSta_phi",
    "GlbMuon_Glb_chi2OverDf","GlbMuon_Tk_chi2OverDf","GlbMuon_Sta_chi2OverDf",
    "TkMuon_chi2OverDf","StaMuon_chi2OverDf",
    "GlbMuon_qComparison",
    "Res_TkGlb_qOverp","Res_GlbSta_qOverp","Res_TkSta_qOverp",
    "Res_TkGlb_oneOverp","Res_GlbSta_oneOverp","Res_TkSta_oneOverp",
    "Res_TkGlb_qOverpt","Res_GlbSta_qOverpt","Res_TkSta_qOverpt",
    "Res_TkGlb_oneOverpt","Res_GlbSta_oneOverpt","Res_TkSta_oneOverpt",
    "StaRh_Frac_inGlb","TkRh_Frac_inGlb",
    "StaRh_inGlb_Div_RhAssoSta","TkRh_inGlb_Div_RhAssoTk",
    "GlbRh_Div_RhAssoStaTk","invalidRh_Frac_inTk"
  };//31

  for(int i = 0 ; i < 5; ++i) {
  e_canvas = newCanvas("mra/"+histoName[i], histoName[i]);
  e_objCol = makeObjectCollection(dirList,histoName[i]);
  e_theLegend = drawObjectCollection(e_objCol,true,legendPt);
  //aaa  e_graphCol = makeTGfromTHCollection(e_objCol);
  //aaa  e_theLegend = drawObjectCollection(e_graphCol,true,legendPt);
  //((TGraph*)e_graphCol->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.9,1.01);
  //((TGraph*)e_graphCol->First())->GetHistogram()->GetYaxis()->SetTitle("Efficiency");
  //((TGraph*)e_graphCol->First())->GetHistogram()->GetXaxis()->SetTitle(e_XAxis[iXaxis]);
  e_theLegend->Draw("same");
  }

  //-------------------------------

  
  //printCanvasesType(".eps");
}
