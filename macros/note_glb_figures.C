TLegend * note_glb_figures()
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

  TList * fileList = makeFileCollection("my2112FileList.txt");

  TString directories[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/MultiTrack/globalMuons_tpToGlbAssociation"
  }

  TList * dirList = makeDirectoryCollection(fileList,directories,1);

  TString directoriesComposites[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/MultiTrack"
  };

  TList * dirListComposites = makeDirectoryCollection(fileList,directoriesComposites,1);

  TString legendPt[] = {"muPt10","muPt100","muPt200","muPt500","muPt1000"}
  //TString legendPt[] = {"muPt100"};

  //----------------------------------------------------
  // Now do each specific figure
  //----------------------------------------------------

  //----------------------------------------------------
  // Algorithm Efficiencies
  //----------------------------------------------------
  
  TString ec_level[] = {"global/glb"};
  TString ec_levelName[] = { "GlobalMuons" };
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
      for (int iXaxis=0;iXaxis < 1; ++iXaxis) {
	e_canvas = newCanvas("FigGLB/efficiencies/"+ec_level[iLevel]+"_"+ec_quantity[iQuantity]+ec_xAxis2[iXaxis],ec_levelName[iLevel]+" "+ec_histo[0]+" "+ec_quantity[iQuantity]+ec_xAxis2[iXaxis]);

	ec_objCol =  makeObjectCollection(dirListComposites,ec_histo[0]+"_"+ec_quantity2[iQuantity]+"_"+ec_xAxis[iXaxis]);
	drawObjectCollection(ec_objCol,false);
	ec_graphCol = makeTGfromTHCollection(ec_objCol);
	drawObjectCollection(ec_graphCol,true,legendPt);
	((TGraph*)ec_graphCol->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.9,1.01);
	((TGraph*)ec_graphCol->First())->GetHistogram()->GetYaxis()->SetTitle("Efficiency");
	((TGraph*)ec_graphCol->First())->GetHistogram()->GetXaxis()->SetTitle(ec_XAxis[iXaxis]);
      }
    }
  }
  

  //----------------------------------------------------
  // Absolute Efficiencies
  //----------------------------------------------------
  
  TString e_level[] = {"global/glb"};
  TString e_levelName[] = { "GlobalMuons" };
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
	e_canvas = newCanvas("FigGLB/efficiencies/"+e_level[iLevel]+"_"+e_quantity[iQuantity]+e_xAxis2[iXaxis],e_levelName[iLevel]+" "+e_histo[0]+" "+e_quantity[iQuantity]+e_xAxis2[iXaxis]);

	e_objCol =  makeObjectCollection(dirList,e_histo[0]+e_xAxis[iXaxis]);
	drawObjectCollection(e_objCol,false);
	e_graphCol = makeTGfromTHCollection(e_objCol);
	drawObjectCollection(e_graphCol,true,legendPt);
	((TGraph*)e_graphCol->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.9,1.01);
	((TGraph*)e_graphCol->First())->GetHistogram()->GetYaxis()->SetTitle("Efficiency");
	((TGraph*)e_graphCol->First())->GetHistogram()->GetXaxis()->SetTitle(e_XAxis[iXaxis]);
      }
    }
  }
  

  //----------------------------------------------------
  // Resolutions
  //----------------------------------------------------
  
  TString level[] = {"global/glb"};
  TString levelName[] = { "GlobalMuons" };
  TString quantity[] = {"res_dxy",
			"res_dz",
			"res_eta",
			"res_phi",
			"res_qpt",
			"res_cotTheta"};
  TString yTitle[] = {"#sigma(#delta d_{xy})[cm]",
		      "#sigma(#delta d_{z}) [cm]",
		      "#sigma(#delta #eta)",
		      "#sigma(#delta #phi)[rad]",
		      "#sigma(#delta (q/p_{t}))",
		      "#sigma(cot #theta)"};
  TString histo[] = {"dxyres",
		     "dzres",
		     "etares",
		     "phires",
		     "ptres",
		     "cotThetares"};
  TString xAxis[]  = {"eta","pt","phi"}; 
  TString XAxis[]  = {"|#eta|","p_{T} (GeV/c)","#phi (rad)"}; 

  int canvasCounter = 01;
  TCanvas * canvas;
  TLegend * theLegend;
  TH1 * h;
  TList * objCol;
  TList * fitCol;
  for (int iLevel=0;iLevel < 1; ++iLevel) {
    for (int iQuantity=0;iQuantity < 6; ++iQuantity) {
      for (int iXaxis=0;iXaxis < 3; ++iXaxis) {
	if(iQuantity==2 && iXaxis==1) continue;
	if(iXaxis==2 && (iQuantity != 3 || iQuantity !=4)) continue;
	canvas = newCanvas("FigParam/resolutions/"+level[iLevel]+"_"+quantity[iQuantity]+"_vs_"+xAxis[iXaxis],"Resolutions "+levelName[iLevel]+" "+quantity[iQuantity]+"_vs_"+xAxis[iXaxis] );

	objCol = makeObjectCollection(dirList,histo[iQuantity]+"_vs_"+xAxis[iXaxis]);
	fitCol = makeFitCollection(objCol,2);
	((TGraph*)fitCol->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.00001,0.5);
	theLegend = drawObjectCollection(fitCol,true,legendPt);
	((TGraph*)fitCol->First())->GetHistogram()->GetYaxis()->SetTitle(yTitle[iQuantity]);
	((TGraph*)fitCol->First())->GetHistogram()->GetXaxis()->SetTitle(XAxis[iXaxis]);
	canvas->SetLogy(1);
	theLegend->SetX1NDC(0.2);
	theLegend->SetX2NDC(0.5);
	theLegend->SetY1NDC(0.7);
	theLegend->SetY2NDC(0.9);
	theLegend->Modify();
	canvas->Update();
      }
    }
  }
  

  //printCanvasesType(".eps");
}
