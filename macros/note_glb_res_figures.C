TLegend * note_glb_res_figures()
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

  TString dirName_("MultiTrack/globalMuons_tpToGlbAssociation");
  TString figDirName_("global/glb");
  TString figLongName_("Global Muons");

  TList * fileList = makeFileCollection("my2112FileList.txt");

  TString directories[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName_
  }

  TList * dirList = makeDirectoryCollection(fileList,directories,1);

  TString legendPt[] = {"muPt10","muPt100","muPt200","muPt500","muPt1000"}
  //TString legendPt[] = {"muPt100"};

  //----------------------------------------------------
  // Now do each specific figure
  //----------------------------------------------------

  //----------------------------------------------------
  // Resolutions
  //----------------------------------------------------
  
  TString level[] = {figDirName_};
  TString levelName[] = { figLongName_ };
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

  printCanvasesType(".eps");
}
