TLegend * note_glb_res_unbinned_figures()
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

  TString dirName_("MultiTrack/general_tpToTkmuAssociation");
  TString dirName_3("MultiTrack/globalMuons_tpToGlbAssociation");
  TString dirName_2("MultiTrack/standAloneMuons_UpdatedAtVtx_tpToStaAssociation");
  TString figDirName_("global/glb");
  TString figLongName_("Muons");

  TList * fileList = makeFileCollection("my2112FileList_unbinned.txt");

  TString directories[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName_,
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName_2,
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName_3
  }

  TList * dirList = makeDirectoryCollection(fileList,directories,3);

  TString legendPt[] = {"General","Stand Alone","Global"}

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
			"res_pt",
			"res_cotTheta",
			"res_p",
			"res_qpt"};
  TString yTitle[] = {"#sigma(#delta d_{xy})[cm]",
		      "#sigma(#delta d_{z}) [cm]",
		      "#sigma(#delta #eta)",
		      "#sigma(#delta #phi)[rad]",
		      "#sigma(#delta p_{T})",
		      "#sigma(cot #theta)",
		      "#sigma(#delta p)",
		      "#sigma(#delta q/p_{T})"};
  TString histo[] = {"dxyres",
		     "dzres",
		     "etares",
		     "phires",
		     "ptres",
		     "cotThetares",
		     "ErrP",
		     "ErrQPt"};
  TString xAxis[]  = {"eta","pt","phi"}; 
  TString xAxis2[]  = {"Eta","Pt","Phi"}; 
  TString XAxis[]  = {"|#eta|","p_{T} (GeV/c)","#phi (rad)"}; 

  int canvasCounter = 01;
  TCanvas * canvas;
  TCanvas * canvas2;
  TLegend * theLegend;
  TH1 * h;
  TList * objCol;
  TList * fitCol;
  for (int iLevel=0;iLevel < 1; ++iLevel) {
    for (int iQuantity=0;iQuantity < 6; ++iQuantity) {
      for (int iXaxis=0;iXaxis < 3; ++iXaxis) {
	if(iXaxis==1 && (iQuantity==2)) continue;
	if(iXaxis==2 && (iQuantity == 0 || iQuantity ==1 || iQuantity == 2 || iQuantity == 5)) continue;

	canvas = newCanvas("FigParam/resolutions/"+quantity[iQuantity]+"_vs_"+xAxis[iXaxis]+"_unbinned","Resolutions "+quantity[iQuantity]+"_vs_"+xAxis[iXaxis] );

	objCol = makeObjectCollection(dirList,histo[iQuantity]+"_vs_"+xAxis[iXaxis]);
	fitCol = makeFitCollection(objCol,2);
	((TGraph*)fitCol->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.00001,1.);
	theLegend = drawObjectCollection(fitCol,true,legendPt);
	((TGraph*)fitCol->First())->GetHistogram()->GetYaxis()->SetTitle(yTitle[iQuantity]);
	((TGraph*)fitCol->First())->GetHistogram()->GetXaxis()->SetTitle(XAxis[iXaxis]);
	canvas->SetLogy(1);
	theLegend->SetX1(0.2);
	theLegend->SetX2(0.5);
	theLegend->SetY1(0.7);
	theLegend->SetY2(0.9);
	theLegend->Modify();
	canvas->Update();
      }
    }
  }

  //------------------------------------
  TString dirName2_("RecoMuon_TrackAssoc/Trk");
  TString dirName2_2("RecoMuon_TrackAssoc/Sta");
  TString dirName2_3("RecoMuon_TrackAssoc/Glb");
  TString directoriesRM[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName2_,
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName2_2,
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName2_3
  }

  TList * dirList2 = makeDirectoryCollection(fileList,directoriesRM,3);

  for (int iLevel=0;iLevel < 1; ++iLevel) {
    for (int iQuantity=6;iQuantity < 8; ++iQuantity) {
      for (int iXaxis=0;iXaxis < 2; ++iXaxis) {
	if(iXaxis>=1 && iQuantity==6) continue;

	canvas = newCanvas("FigParam/resolutions/"+level[iLevel]+"_"+quantity[iQuantity]+"_vs_"+xAxis[iXaxis],"Resolutions "+levelName[iLevel]+" "+quantity[iQuantity]+"_vs_"+xAxis[iXaxis]);
	
	objCol = makeObjectCollection(dirList2,histo[iQuantity]+"_vs_"+xAxis2[iXaxis]);

	fitCol = makeFitCollection(objCol,2);
	((TGraph*)fitCol->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.00001,1.);
	((TGraph*)fitCol->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.,500.);
	theLegend = drawObjectCollection(fitCol,true,legendPt);
	((TGraph*)fitCol->First())->GetHistogram()->GetYaxis()->SetTitle(yTitle[iQuantity]);
	((TGraph*)fitCol->First())->GetHistogram()->GetXaxis()->SetTitle(XAxis[iXaxis]);
	canvas->SetLogy(1);
	theLegend->SetX1(0.2);
	theLegend->SetX2(0.5);
	theLegend->SetY1(0.7);
	theLegend->SetY2(0.9);
	theLegend->Modify();
	canvas->Update();
      }
    }
  }
  
  //printCanvasesType(".eps");
}
