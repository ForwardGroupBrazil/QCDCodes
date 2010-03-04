TLegend * note_mtv_figs()
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

  //TString dirName_("MultiTrack/globalMuons_tpToGlbAssociation");
  //TString level = "glb";
  TString dirName_("MultiTrack/standAloneMuons_UpdatedAtVtx_tpToStaUpdAssociation");
  TString level = "staUpd";

  TList * fileList = makeFileCollection("MyBinnedList.txt");

  TString directories[] = {
    "/DQMData/Run 1/Muons/Run summary/RecoMuonV/"+dirName_
  }

  TList * dirList = makeDirectoryCollection(fileList,directories,1);


  TString legendPt[] = {"muPt10","muPt100","muPt500","muPt1000","muPt0_500"}

  //----------------------------------------------------
  // Now do each specific figure
  //----------------------------------------------------

  TString histoName[] = {"effic","efficPt","effic_vs_phi","effic_vs_hit",
			 //"nhits_vs_eta_pfx",
			 "hits","chi2","chi2_prob",
			 //"chi2_vs_eta_pfx",
			 "pullPt","pullQoverp",
			 "pullPhi","pullTheta",
			 "pullDxy","pullDz",
			 //
			 //"ptres_vs_eta","etares_vs_eta",
			 //"phires_vs_eta","cotThetares_vs_eta",
			 //"dxyres_vs_eta","dzres_vs_eta"
			 //
			 "phires_vs_eta_Sigma","cotThetares_vs_eta_Sigma",
			 "dxyres_vs_eta_Sigma","dzres_vs_eta_Sigma",
			 "ptres_vs_eta_Sigma",
			 "phires_vs_eta_Mean","cotThetares_vs_eta_Mean",
			 "dxyres_vs_eta_Mean","dzres_vs_eta_Mean",
			 "ptres_vs_eta_Mean",
			 "phires_vs_pt_Sigma","cotThetares_vs_pt_Sigma",
			 "dxyres_vs_pt_Sigma","dzres_vs_pt_Sigma",
			 "ptres_vs_pt_Sigma"

  };//28

  TString xTitle[] = {"#eta","p_{T} [GeV]","#phi [rad]","nHits",
		      "nHits","#chi^{2}","#chi^{2}/DoF",
		      "","",
		      "","",
		      "","",
		      "#eta","#eta",
		      "#eta","#eta",
		      "#eta",
		      "#eta","#eta",
		      "#eta","#eta",
		      "#eta",
		      "p_{T} [GeV]","p_{T} [GeV]",
		      "p_{T} [GeV]","p_{T} [GeV]",
		      "p_{T} [GeV]",
  };

  TString yTitle[] = {"Efficiency","Efficiency","Efficiency","Efficiency",
		      "Entries","Entries","Entries",
		      "Entries","Entries",
		      "Entries","Entries",
		      "Entries","Entries",
		      "#sigma(#delta #phi) [rad]","#sigma(#delta cot(#theta))",
		      "#sigma(#delta d_{xy}) [cm]","#sigma(#delta z_{0}) [cm]",
		      "#sigma(#delta p_{t}/p_{t}) ",
		      "#delta #phi [rad]","#delta cot(#theta)",
		      "#delta d_{xy} [cm]","#delta z_{0} [cm]",
		      "#delta p_{T}/p_{T}",
		      "#sigma(#delta #phi) [rad]","#sigma(#delta cot(#theta))",
		      "#sigma(#delta d_{xy}) [cm]","#sigma(#delta z_{0}) [cm]",
		      "#sigma(#delta p_{t}/p_{t}) ",
  };

  for(int i =0; i < 28; ++i) {
    e_canvas = newCanvas("mtv/"+level+"/"+histoName[i], level + " " + histoName[i],0,0,2);
    e_objCol = makeObjectCollection(dirList,histoName[i]);
    e_theLegend = drawObjectCollection(e_objCol,true,legendPt);
    //aaa  e_graphCol = makeTGfromTHCollection(e_objCol);
    //aaa  e_theLegend = drawObjectCollection(e_graphCol,true,legendPt);
    //((TGraph*)e_graphCol->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.9,1.01);
    ((TH1*)e_objCol->First())->GetYaxis()->SetTitle(yTitle[i]);
    ((TH1*)e_objCol->First())->GetXaxis()->SetTitle(xTitle[i]);
    e_theLegend->Draw("same");
  }

  //-------------------------------

  
  //printCanvasesType(".eps");
}
