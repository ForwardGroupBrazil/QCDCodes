void DrawHistograms(TString HISTO,TString XTITLE,float YMAX,float XMIN,float XMAX,int REBIN,bool LOGY,bool PRINT)
{
  gROOT->ForceStyle();
  const int N = 11;
  TFile *inf[N];
  TH1F *h[N];
  TString PATH("/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/");
  TString fileName[N] = {
    "flatTree_data_preselect_hard_tmva_histos.root",
    "flatTree_VBF-Powheg125_preselect_hard_tmva_histos.root",
    "flatTree_GluGlu-Powheg125_preselect_hard_tmva_histos.root",
    "flatTree_T_preselect_hard_tmva_histos.root",
    "flatTree_Tbar_preselect_hard_tmva_histos.root",
    "flatTree_TTJets_preselect_hard_tmva_histos.root",
    "flatTree_ZJets_preselect_hard_tmva_histos.root",
    "flatTree_QCD-HT100_preselect_hard_tmva_histos.root",
    "flatTree_QCD-HT250_preselect_hard_tmva_histos.root",
    "flatTree_QCD-HT500_preselect_hard_tmva_histos.root",
    "flatTree_QCD-HT1000_preselect_hard_tmva_histos.root"
  };
  float LUMI = 18530;
  float XSEC[N]    = {1,0.911,11.26,71.29,43.56,225.197,650,1.036e+7,2.76e+5,8426,204.};
  int   NEVENTS[N] = {1,991754,999908,3748227,1935072,6923750,8130192,50023348,26981684,30522164,13808234};
  int FILLCOLOR[N] = {kBlack,kOrange,kOrange+4,kGreen+2,kGreen+3,kGreen-8,kMagenta,kBlue-6,kBlue-8,kBlue-9,kBlue-10};
  int LINECOLOR[N] = {kBlack,kOrange,kOrange+4,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack};
  int LINESIZE[N]  = {1,2,2,1,1,1,1,1,1,1};
  THStack *hs = new THStack("bkg","bkg");
  TH1F *hBkg;
  TH1F *hQCD;
  TH1F *hTop;
  for(int iFile=0;iFile<N;iFile++) {
    inf[iFile] = TFile::Open(PATH+fileName[iFile]);
    h[iFile] = (TH1F*)inf[iFile]->Get(HISTO);
    h[iFile]->Sumw2();
    h[iFile]->Rebin(REBIN);
    h[iFile]->SetLineColor(LINECOLOR[iFile]);
    h[iFile]->SetLineWidth(LINESIZE[iFile]);
    if (iFile == 0) {
      hBkg = (TH1F*)h[iFile]->Clone("hBkg");
      hBkg->Reset();
      hQCD = (TH1F*)h[iFile]->Clone("hQCD");
      hQCD->Reset();
      hTop = (TH1F*)h[iFile]->Clone("hTop");
      hTop->Reset(); 
    }
    else {
      float wt = LUMI*XSEC[iFile]/NEVENTS[iFile];
      h[iFile]->Scale(wt);
      if (iFile > 2) {
        h[iFile]->SetFillColor(FILLCOLOR[iFile]);
        hBkg->Add(h[iFile]);
      }  
      if (iFile > 2 && iFile < 5) {
        hTop->Add(h[iFile]);
      } 
      if (iFile > 6) {
        hQCD->Add(h[iFile]);
      }
    }
  }  
  float NQCD = hQCD->Integral();
  float NTOT = hBkg->Integral();
  float NDAT = h[0]->Integral();
  float kfactor = (NDAT - NTOT)/NQCD + 1;
  cout<<"k-factor: "<<kfactor<<endl;
  /*
  for(int iFile=6;iFile<N;iFile++) {
    if (iFile > 6) {
      h[iFile]->Scale(kfactor);
    }  
  }
  */
  hQCD->Scale(kfactor);
  hBkg->Scale(1+(kfactor-1)*NQCD/NTOT);
  TH1F *hRatio = (TH1F*)h[0]->Clone("Ratio");
  hRatio->Divide(hBkg);
  
  hTop->SetFillColor(kGreen+2);
  hQCD->SetFillColor(kBlue-8);

  hs->Add(hTop);
  hs->Add(h[5]);
  hs->Add(h[6]);
  hs->Add(hQCD);

  TCanvas *can = new TCanvas("can_"+HISTO,"can_"+HISTO,900,600);
  can->cd(1);
  can->SetBottomMargin(0.3);
  can->SetRightMargin(0.2);
  if (LOGY) gPad->SetLogy();
  TH1F *haux = (TH1F*)hBkg->Clone("aux");
  haux->Reset();
  haux->GetYaxis()->SetTitle("Events");
  
  //haux->GetXaxis()->SetTitle(XTITLE);
  haux->GetXaxis()->SetLabelSize(0.0);
  haux->GetYaxis()->SetNdivisions(505);
  haux->GetXaxis()->SetNdivisions(505);
  haux->GetXaxis()->SetRangeUser(XMIN,XMAX);
  haux->SetMinimum(0.5);
  //haux->SetMaximum(1.2*h[0]->GetBinContent(h[0]->GetMaximumBin()));
  haux->SetMaximum(YMAX);
  haux->Draw();
  hs->Draw("same hist");
  h[0]->Draw("same E");
  h[1]->Draw("same hist");
  h[2]->Draw("same hist");
  TLegend *leg = new TLegend(0.8,0.5,0.92,0.9);
  leg->SetHeader("Presel. & Trigger");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(h[0],"Data (18.5 fb^{-1})","P");
  leg->AddEntry(hQCD,"QCD","F");
  leg->AddEntry(h[6],"ZJets","F");
  leg->AddEntry(h[5],"TTJets","F");
  leg->AddEntry(hTop,"Single Top","F");
  leg->AddEntry(h[1],"Signal VBF","L");
  leg->AddEntry(h[2],"Signal GG","L"); 
  leg->Draw();
  gPad->RedrawAxis();
  TPad* pad = new TPad("pad", "pad", 0., 0., 1., 1.);
  pad->SetTopMargin(0.7);
  pad->SetRightMargin(0.2);
  pad->SetFillColor(0);
  pad->SetFillStyle(0);
  pad->Draw();
  pad->cd(0);
  gPad->SetGridy();
  hRatio->SetMarkerSize(0.7);
  hRatio->SetMinimum(0);
  hRatio->SetMaximum(2);
  hRatio->GetYaxis()->SetTitle("Data / MC");
  hRatio->GetXaxis()->SetTitle(XTITLE);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetTickLength(0.06);
  hRatio->GetYaxis()->SetTitleSize(0.04);
  hRatio->GetYaxis()->SetLabelSize(0.03);
  hRatio->GetXaxis()->SetNdivisions(505);
  hRatio->GetXaxis()->SetRangeUser(XMIN,XMAX);
  hRatio->Draw();
  
  if (PRINT) {
    can->Print(TString(can->GetTitle())+".gif");
  }
  
}















