//
//
int gLineColor=1;
TLegend *gLegend;
int gOpt=1;
int gLineStyle=1;

bool gForce= false;


TH1* chargeId(TString sampleLoc, TString recoStep, TString histoType,int misOp = 0){
  //gROOT->LoadMacro("macros.C");     // Load service macros

  gLineColor=1;
  gLineStyle=1;
  
  gLegend = new TLegend(.40, .15, 0.63, .35, ""); // eff 5 GeV seed

  //  TString recoStep = "StandAloneMuonsUpdatedAtVtx";
  TH1* hhh = getInclusive(sampleLoc,recoStep,histoType,misOp);

  hhh->SetMinimum(0.9);
  hhh->SetMaximum(1.01);

  gLegend->SetFillColor(0);
  gLegend->Draw("same");
  
  return hhh;

}
  

TH1 *get(TString sampleLoc, TString sample, TString recoStep, TString histoType, int misOp = 0){
  
  TString histoName = recoStep+"ResoCharge"+histoType;
  
  TH2F *hSta = (TH2F*)getHisto(sampleLoc,sample,histoName);

  int y1 = 2;
  int y2 = -1;
  
  if(misOp==1){
    int y1 = 0;
    int y2 = 1;
  }
  // 2,-1 --> correct ass prob
  // 0,1 --> mis-ass prob
  TH1D *hStaPr = hSta->ProjectionX("_",y1,y2);
  
  TString histoTypeSim = histoType;
  histoTypeSim.ReplaceAll("Vs","");
  TH1F *hSim = (TH1F*)getHisto(sampleLoc,sample,recoStep+"SimAsso"+histoTypeSim);
  TH1F *hStaPr2 =  hStaPr->Clone();
  //  hStaPr2->Sumw2();
  // hStaPr2->Divide(hSim);
  TH1 *hres = computeEfficiency(hStaPr2,hSim);  
  drawHisto(hres);
  TString type = reformatLegend(sample);
  gLegend->AddEntry(hres, type, "PL");

  return hres;
}


TH1F *get(TString sampleLoc, TString sample, TString recoStep, TString histoType, int misOp,
	  int msty, int lsty){

  TH1F *histo=get(sampleLoc,sample,recoStep,histoType,misOp);
  histo->SetLineStyle(lsty);
  histo->SetMarkerStyle(msty);
  gLineStyle = lsty;
  TH1F *histo2=histo->Clone();
  histo2->SetLineStyle(1);
  //  --gLineColor;
  if(--gLineColor==5) --gLineColor;
  drawHisto(histo2);
  histo2->SetMarkerSize(msty ==21 ? 0.7 : 0.1);
  return histo2;
}



TH1F *getInclusive(TString &sampleLoc, TString recoStep, TString &histoType,int misOp = 0){
  //  gLineColor = 2;
  // gLineStyle = 10;
  TH1F *histo = get(sampleLoc,"mu10",recoStep,histoType,misOp,21,1);
  histo->SetMarkerSize(0.7);
  
  //histo->GetXaxis()->SetTitle("#eta");
  //histo->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //histo->GetXaxis()->SetTitle("#phi (rad)");
  //  histo->GetXaxis()->SetTitleSize(0.065); // eta
  // histo->GetXaxis()->SetTitleOffset(0.55); //eta

  //histo->GetXaxis()->SetTitleSize(0.05); // phi
  //histo->GetXaxis()->SetTitleOffset(0.75); //phi
  
  histo->GetYaxis()->SetTitle("Charge id probability");
  //histo->GetYaxis()->SetTitleOffset(1.37);
  //histo->GetYaxis()->SetNdivisions(511);
  //histo->GetYaxis()->SetLabelSize(0.04);

  //histo->Clear();




  //  gLineColor = 1;
  // gLineStyle = 1;
  // change the requirement in macro.C in order to have "same" option in Draw()

  //  get(sampleLoc,"mu10",histoName,21,2);
  //get(sampleLoc,"mu50",recoStep,histoType,misOp,22,3);
  get(sampleLoc,"mu100",recoStep,histoType,misOp,23,4);
  get(sampleLoc,"mu200",recoStep,histoType,misOp,24,5);
  get(sampleLoc,"mu500",recoStep,histoType,misOp,25,6);
  get(sampleLoc,"mu1000",recoStep,histoType,misOp,26,7);
  return histo;
}

TH1F *computeEfficiency(const TH1F *num, const TH1F *denom){
  
  //  num->Rebin(5); denom->Rebin(5);

  TH1F *hEff = (TH1F*) num->Clone();
  
  hEff->Divide(denom);

  // Set the error accordingly to binomial statistics
  int nBinsEta = hEff->GetNbinsX();
  for(int bin = 1; bin <=  nBinsEta; bin++) {
    float nDenomHit = denom->GetBinContent(bin);
    float eff = hEff->GetBinContent(bin);
    float error = 0;
    if(nDenomHit != 0 && eff <= 1) {
      error = sqrt(eff*(1-eff)/nDenomHit);
    }
    hEff->SetBinError(bin, error);
  }
  return hEff;
}
