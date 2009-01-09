TLegend * note_trk_pull_unbinned_figures()
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
  TString figDirName_("tracker/tk");
  TString figLongName_("General Tracks");

  TList * fileList = makeFileCollection("my2112FileList_unbinned.txt");

  TString directories[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName_
  }

  TList * dirList = makeDirectoryCollection(fileList,directories,1);

  TString legendPt[] = {"muPt0_500"}


  //----------------------------------------------------
  // Now do each specific figure
  //----------------------------------------------------


  TString p_level[] = {figDirName_};
  TString p_levelName[] = {figLongName_};
  
  TString p_quantity[] = {"Pull_dxy",
			  "Pull_dz",
			  "Pull_theta",
			  "Pull_phi",
			  "Pull_qp",
			  "Pull_pt"
  };
  TString p_yTitle[] = {"#sigma_{pull}(d_{xy})",
			"#sigma_{pull}(d_{z})",
			"#sigma_{pull}(#theta)",
			"#sigma_{pull}(#phi)",
			"#sigma_{pull}(q/p)",
			"#sigma_{pull}(p_{t})"
  };
  TString p_histo[] = {"pullDxy",
		       "pullDz",
		       "pullTheta",
		       "pullPhi",
		       "pullQoverp",
		       "pullPt"
  };
  //  TString p_xAxis[]  = {"eta","phi"}; 
  //  TString p_XAxis[]  = {"Eta","Phi"}; 

  TCanvas * p_canvas;
  TLegend * p_theLegend;
  TH1 * p_h;
  TList * p_objCol;
  TList * p_fitCol;
  for (int iLevel=0;iLevel < 1; ++iLevel) {
    for (int iQuantity=0; iQuantity < 6; ++iQuantity) {
      //for (int iXaxis=0;iXaxis < 2; ++iXaxis) {
      //if(iXaxis==1 && (iQuantity==0 || iQuantity==1)) continue;

	p_canvas = newCanvas("FigParam/pulls/"+p_level[iLevel]+"_"+p_quantity[iQuantity]+"_unbinned","Pulls "+p_levelName[iLevel]+" "+p_quantity[iQuantity] );
	p_objCol = makeObjectCollection(dirList,p_histo[iQuantity]);
	makeFitCollection(p_objCol,2);
	//	((TGraph*)p_fitCol->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.,3.);
	//	p_theLegend = drawObjectCollection(p_fitCol,true);
	p_theLegend = drawObjectCollection(p_objCol,true);
	((TH1*)p_objCol->First())->GetYaxis()->SetTitle(p_yTitle[iQuantity]);
	//	((TH1*)p_objCol->First())->GetXaxis()->SetTitle(p_XAxis[iXaxis]);
	p_canvas->SetLogy(0);
	p_theLegend->SetX1NDC(0.2);
	p_theLegend->SetX2NDC(0.5);
	p_theLegend->SetY1NDC(0.7);
	p_theLegend->SetY2NDC(0.9);
	p_theLegend->Modify();
	p_theLegend->Delete();
	p_canvas->Update();
	//}
    }
  }

    //printCanvasesType(".eps");
}
