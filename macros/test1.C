void test1(TString &dir="./"){
  gROOT->LoadMacro("adamStyles.C");

  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();

  gROOT->LoadMacro("adamCanvasMacros.C");
  gROOT->LoadMacro("adamDrawMacros_old.C");
  gROOT->LoadMacro("adamDrawLoop.C");

  gROOT->LoadMacro("adamGetObjMacros.C");
  gROOT->LoadMacro("adamDrawMacros.C");

  TSortedList * directories = new TSortedList();
  directories->Add(getDirectory("DQM_file.root","/DQMData/Run 1/RecoMuonV/Run summary/MultiTrack/general_tpToTkmuAssociation")   );
  directories->Add(getDirectory("DQM_file.root","/DQMData/Run 1/RecoMuonV/Run summary/MultiTrack/globalMuons_tpToGlbAssociation")   );

  TString histos[] = {
    "efficPt", "effic"
  }


  //  TIter iterH(histos);

  TDirectory* c;
  TString* h;
  TH1 * tmp;
  TGraph * tmpG;


  for (int i = 0; i<2 ; ++i) {
    cout << i << endl;
  TCanvas * c1_n1 = new TCanvas(histos[i],histos[i]);
  TCanvas * c1_n2 = new TCanvas(histos[i]+"graph",histos[i]+"graph");

  TIter iter(directories);
  int color = 1;
  while( (c = (TDirectory *)iter()) )
    {
      cout << histos[i] << " " << c->GetName() << endl;
      TString optColorH(color==1 ? "" : "same");
      TString optErrH = ("E1X0") + optColorH;

      TString optColorG(color==1 ? "APL" : "PL same");
      cout << color << " " << optColorG.Data() << endl;
      TString optErrG = (" ")+optColorG;
      cout << color << " " << optErrG.Data() << endl;

      //tmp = getHistogram(c,"effic");
      tmp = getHistogram(c,histos[i]);
      tmp->SetMarkerColor(color);
      tmp->SetLineColor(color);

      c1_n1->cd();
      tmp->Draw(optErrH);

      c1_n2->cd();
      tmpG = convertTH1toTGraph(tmp);
      tmpG->Draw(optErrG);

      color++;
    }

  }

}
