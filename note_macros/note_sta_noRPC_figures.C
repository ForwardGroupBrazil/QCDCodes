TLegend * note_sta_noRPC_figures()
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

  TString dirName_("MultiTrack/standAloneMuons_UpdatedAtVtx_tpToStaAssociation");
  TString dirName_2("MultiTrack/standAloneMuonsNoRPC_UpdatedAtVtx_tpToStaMuonNoUpAssociation");
  TString figDirName_("standAlone/sta");
  TString figLongName_("Stand Alone Muons");
  TString figDirName2_("FigSTA");

  bool doAlgEff_ = false;

  TList * fileList = makeFileCollection("my2112FileList_unbinned.txt");

  TString directories[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName_,
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName_2
  }

  TList * dirList = makeDirectoryCollection(fileList,directories,2);

  TString directoriesComposites[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/MultiTrack"
  };

  TList * dirListComposites = makeDirectoryCollection(fileList,directoriesComposites,1);

  TString legendPt[] = {"Stand Alone","Stand Alone: no RPC"}

  //----------------------------------------------------
  // Now do each specific figure
  //----------------------------------------------------

  //----------------------------------------------------
  // Algorithm Efficiencies
  //----------------------------------------------------
  if(doAlgEff_) {  
  TString ec_level[] = {figDirName_};
  TString ec_levelName[] = { figLongName_ };
  TString ec_quantity[] = {"sta_eff","tk_eff"};
  TString ec_quantity2[] = {"GlbSta","GlbTk"};
  TString ec_histo[] = {"Eff"};
  TString ec_xAxis[] = {"Eta","Hit","Pt"};
  TString ec_xAxis2[] = {"_vs_eta","_vs_hit","_vs_pt"};
  TString ec_XAxis[] = {"|#eta|","n Hits","p_{T} (GeV)"};

  int ec_canvasCounter = 01;
  TCanvas * ec_canvas;
  TLegend * ec_theLegend;
  TH1 * ec_h;
  TList * ec_objCol;
  TList * ec_graphCol;

  for (int iLevel=0;iLevel < 1; ++iLevel) {
    for (int iQuantity=0;iQuantity < 2; ++iQuantity) {
      for (int iXaxis=0;iXaxis < 3; ++iXaxis) {
	e_canvas = newCanvas(figDirName2_+"/efficiencies/"+ec_level[iLevel]+"_"+ec_quantity[iQuantity]+ec_xAxis2[iXaxis],ec_levelName[iLevel]+" "+ec_histo[0]+" "+ec_quantity[iQuantity]+ec_xAxis2[iXaxis]);

	ec_objCol =  makeObjectCollection(dirListComposites,ec_histo[0]+"_"+ec_quantity2[iQuantity]+"_"+ec_xAxis[iXaxis]);
	drawObjectCollection(ec_objCol,false);
	ec_graphCol = makeTGfromTHCollection(ec_objCol);
	drawObjectCollection(ec_graphCol,true,legendPt);
	((TGraph*)ec_graphCol->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.5,1.05);
	((TGraph*)ec_graphCol->First())->GetHistogram()->GetYaxis()->SetTitle("Efficiency");
	((TGraph*)ec_graphCol->First())->GetHistogram()->GetXaxis()->SetTitle(ec_XAxis[iXaxis]);
      }
    }
  }
  }  

  //----------------------------------------------------
  // Absolute Efficiencies
  //----------------------------------------------------
  
  TString e_level[] = {figDirName_};
  TString e_levelName[] = { figLongName_ };
  TString e_quantity[] = {"sim_eff"};
  TString e_histo[] = {"effic"};
  TString e_xAxis[] = {"","Pt","_vs_hit"};
  TString e_xAxis2[] = {"_vs_eta","_vs_pt","_vs_hit"};
  TString e_XAxis[] = {"|#eta|","p_{T} (GeV)","n Hits"};

  int e_canvasCounter = 01;
  TCanvas * e_canvas;
  TLegend * e_theLegend;
  TH1 * e_h;
  TList * e_objCol;
  TList * e_graphCol;

  for (int iLevel=0;iLevel < 1; ++iLevel) {
    for (int iQuantity=0;iQuantity < 1; ++iQuantity) {
      for (int iXaxis=0;iXaxis < 3; ++iXaxis) {
	e_canvas = newCanvas(figDirName2_+"/efficiencies/"+e_level[iLevel]+"_"+e_quantity[iQuantity]+e_xAxis2[iXaxis]+"_rpc",e_levelName[iLevel]+" "+e_histo[0]+" "+e_quantity[iQuantity]+e_xAxis2[iXaxis]);

	e_objCol =  makeObjectCollection(dirList,e_histo[0]+e_xAxis[iXaxis]);
	drawObjectCollection(e_objCol,false);
	e_graphCol = makeTGfromTHCollection(e_objCol);
	drawObjectCollection(e_graphCol,true,legendPt);
	((TGraph*)e_graphCol->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.5,1.05);
	((TGraph*)e_graphCol->First())->GetHistogram()->GetYaxis()->SetTitle("Efficiency");
	((TGraph*)e_graphCol->First())->GetHistogram()->GetXaxis()->SetTitle(e_XAxis[iXaxis]);
      }
    }
  }

  //-------------------------------
  TString dirNameRM_("RecoMuon_TrackAssoc/Sta");
  TString dirNameRM_2("RecoMuon_NoRPCTrackAssoc/Sta");
  TString directoriesRM[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirNameRM_,
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirNameRM_2
  }

  TList * dirListRM = makeDirectoryCollection(fileList,directoriesRM,2);

  TString erm_level[] = {figDirName_};
  TString erm_levelName[] = { figLongName_ };
  TString erm_quantity[] = {"sim_eff","sim_eff","sim_eff"};
  TString erm_histo[] = {"EffPhi","MisQProbPt","MisQProbEta"};
  TString erm_xAxis[] = {"phi","pt","eta"};
  TString erm_xAxis2[] = {"_vs_phi","_vs_misQPt","_vs_misQEta"};
  TString erm_XAxis[] = {"#phi [rad]","p_{T} (GeV)","|#eta|"};

  int erm_canvasCounter = 01;
  TCanvas * erm_canvas;
  TLegend * erm_theLegend;
  TH1 * erm_h;
  TList * erm_objCol;
  TList * erm_graphCol;

  for (int iLevel=0;iLevel < 1; ++iLevel) {
    for (int iQuantity=0;iQuantity < 3; ++iQuantity) {
      //      for (int iXaxis=0;iXaxis < 2; ++iXaxis) {

	erm_canvas = newCanvas(figDirName2_+"/efficiencies/"+erm_level[iLevel]+"_"+erm_quantity[iQuantity]+erm_xAxis2[iQuantity]+"_rpc",erm_levelName[iLevel]+" "+erm_histo[iQuantity]+" "+erm_quantity[iQuantity]+erm_xAxis2[iQuantity]);

	erm_objCol =  makeObjectCollection(dirListRM,erm_histo[iQuantity]);
	drawObjectCollection(erm_objCol,false);
	erm_graphCol = makeTGfromTHCollection(erm_objCol);
	drawObjectCollection(erm_graphCol,true,legendPt);
	((TGraph*)erm_graphCol->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.5,1.05);
	((TGraph*)erm_graphCol->First())->GetHistogram()->GetYaxis()->SetTitle("Efficiency");
	((TGraph*)erm_graphCol->First())->GetHistogram()->GetXaxis()->SetTitle(erm_XAxis[iQuantity]);

	//      }
    }
  }

  //-------------------------------

  
  //printCanvasesType(".eps");
}
