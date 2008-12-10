TLegend * note_glb_pull_figures()
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


  TString p_level[] = {figDirName_};
  TString p_levelName[] = {figLongName_};
  
  TString p_quantity[] = {"pull_dxy",
			"pull_dz",
			"pull_theta",
			"pull_phi",
			"pull_qpt"
			
			};
  TString p_yTitle[] = {"#sigma_{pull}(#delta d_{xy})",
		      "#sigma_{pull}(#delta d_{z})",
		      "#sigma_{pull}(#delta #theta)",
		      "#sigma_{pull}(#delta #phi)",
		      "#sigma_{pull}((q/p_{t}))"
		      
		      };
  TString p_histo[] = {"dxypull",
		     "dzpull",
		     "thetapull",
		     "phipull",
		     "ptpull"
		     
		     };
  TString p_xAxis[]  = {"eta"}; 
  TString p_XAxis[]  = {"Eta"}; 

  TCanvas * p_canvas;
  TLegend * p_theLegend;
  TH1 * p_h;
  TList * p_objCol;
  TList * p_fitCol;
  for (int iLevel=0;iLevel < 1; ++iLevel) {
    for (int iQuantity=0;iQuantity < 5; ++iQuantity) {
      for (int iXaxis=0;iXaxis < 1; ++iXaxis) {
	p_canvas = newCanvas("FigParam/pulls/"+p_level[iLevel]+"_"+p_quantity[iQuantity]+"_vs_"+p_xAxis[iXaxis],"Pulls "+p_levelName[iLevel]+" "+p_quantity[iQuantity]+"_vs_"+p_xAxis[iXaxis] );
	p_objCol = makeObjectCollection(dirList,p_histo[iQuantity]+"_vs_"+p_xAxis[iXaxis]);
	p_fitCol = makeFitCollection(p_objCol,2);
	((TGraph*)p_fitCol->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.,3.);
	p_theLegend = drawObjectCollection(p_fitCol,true,legendPt);
	((TGraph*)p_fitCol->First())->GetHistogram()->GetYaxis()->SetTitle(p_yTitle[iQuantity]);
	((TGraph*)p_fitCol->First())->GetHistogram()->GetXaxis()->SetTitle(p_XAxis[iXaxis]);
	p_canvas->SetLogy(0);
	p_theLegend->SetX1NDC(0.2);
	p_theLegend->SetX2NDC(0.5);
	p_theLegend->SetY1NDC(0.7);
	p_theLegend->SetY2NDC(0.9);
	p_theLegend->Modify();
	p_canvas->Update();
      }
    }
  }

  printCanvasesType(".eps");
}
