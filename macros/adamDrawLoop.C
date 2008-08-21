//
TLegend *gLegend;

int gLineColor;
int gLineStyle;

int gOpt;
int gType;

bool gForce= false;

/*
 * opt = 1: draw option "E1XO L"
 *
 * type = 1: efficiency
 * type = 2: resolution
 * type = 3: pull
 */

TH1F *adamDrawLoop(TString & sampleLoc,TString & histoName, int iSwitch =0, int opt = 0, int type = 1){
  gLineColor=1;
  gLineStyle=1;

  gOpt = opt;
  gType = type;

  gLegend = new TLegend(.57, .15, 0.74, .32, ""); // eff 5 GeV seed
  //  gLegend = new TLegend(.60, .25, 0.83, .45, ""); // eff 5 GeV seed
  // gLegend = new TLegend(.65, .25, 0.82, .42, "");
  //  gLegend = new TLegend(.60, .25, 0.90, .60, ""); // eff 5 GeV seed

  gLegend->SetFillColor(0);
  
  TH1F *histo = 0;

  switch(iSwitch){
  case 0:{
    histo = getInclusive(sampleLoc,histoName);
    break;
  }
    
  case 1:{
    histo = getUnbinned(sampleLoc,histoName);
    break;
  }

  case 2:{
    histo = getUnbinnedAll(sampleLoc,histoName);
    break;
  }

  case 5:{
    histo = getAndDrawSingleHisto(sampleLoc,"mu5",histoName);
    break;
  }
  case 10:{
    histo = getAndDrawSingleHisto(sampleLoc,"mu10",histoName);
    break;
  }
  case 50:{
    histo = getAndDrawSingleHisto(sampleLoc,"mu50",histoName);
    break;
  }
  case 100:{
    histo = getAndDrawSingleHisto(sampleLoc,"mu100",histoName);
    break;
  }
  case 200:{
    histo = getAndDrawSingleHisto(sampleLoc,"mu200",histoName);
    break;
  }
  case 500:{
    histo = getAndDrawSingleHisto(sampleLoc,"mu500",histoName);
    break;
  }
  case 1000:{
    histo = getAndDrawSingleHisto(sampleLoc,"mu1000",histoName);
    break;
  }
  case 2000:{
    histo = getAndDrawSingleHisto(sampleLoc,"mu2000",histoName);
    break;
  }
  case 3000:{
    histo = getAndDrawSingleHisto(sampleLoc,"mu3000",histoName);
    break;
  }
  case 999:{
    histo = getAndDrawSingleHisto(sampleLoc,"mu0_500",histoName);
    break;
  }
    
  default: {
    histo = getAndDrawSingleHisto(sampleLoc,"mu5",histoName);
    histo = getAndDrawSingleHisto(sampleLoc,"mu10",histoName);
    histo = getAndDrawSingleHisto(sampleLoc,"mu50",histoName);
    histo = getAndDrawSingleHisto(sampleLoc,"mu100",histoName);
    histo = getAndDrawSingleHisto(sampleLoc,"mu200",histoName);
    histo = getAndDrawSingleHisto(sampleLoc,"mu500",histoName);
    histo = getAndDrawSingleHisto(sampleLoc,"mu1000",histoName);
    break;
  }}

  histo->GetXaxis()->SetTitle("#eta");
  //  histo->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //histo->GetXaxis()->SetTitle("#phi (rad)");
  histo->GetXaxis()->SetTitleSize(0.055); // eta
  histo->GetXaxis()->SetTitleOffset(0.55); //eta

  //histo->GetXaxis()->SetTitleSize(0.05); // phi
  //histo->GetXaxis()->SetTitleOffset(0.75); //phi
  


  histo->GetYaxis()->SetTitle("Efficiency");
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->GetYaxis()->SetNdivisions(511);
  histo->GetYaxis()->SetLabelSize(0.04);

  histo->Clear();

  gLegend->Draw("same");
  return histo;
}

TH1F *getInclusive(TString &sampleLoc, TString &histoName){
  //  TH1F *h5 = getAndDrawSingleHisto(sampleLoc,"mu5",histoName);
  //  h5->SetMarkerSize(0.7);

  //  TH1F *h10 = getAndDrawSingleHisto(sampleLoc,"mu10",histoName);
  //  reformatHisto(h10,21,2);

  TH1F *h100 = (gType==1) ? getAndDrawSingleHisto(sampleLoc,"mu100",histoName) : getAndDrawSlicedHisto(sampleLoc,"mu100",histoName,2);
  reformatHisto(h100,22,3);

  //  TH1F *h200 = getAndDrawSingleHisto(sampleLoc,"mu200",histoName);
  //  reformatHisto(h200,23,4);

  //  TH1F *h500 = getAndDrawSingleHisto(sampleLoc,"mu500",histoName);
  //  reformatHisto(h500,24,5);

  TH1F *h1000 = (gType==1) ? getAndDrawSingleHisto(sampleLoc,"mu1000",histoName) : getAndDrawSlicedHisto(sampleLoc,"mu1000",histoName,2);
  reformatHisto(h1000,25,6);

  //  TH1F *h2000 = getAndDrawSingleHisto(sampleLoc,"mu2000",histoName);
  //  reformatHisto(h2000,26,7);

  //  TH1F *h3000 = getAndDrawSingleHisto(sampleLoc,"mu3000",histoName);
  //  reformatHisto(h3000,27,8);

  return h100;
}


TH1F *getUnbinned(TString &sampleLoc, TString &histoName){
  TH1F *h = getAndDrawSingleHisto(sampleLoc,"mu0_500",histoName);
  h->SetMarkerSize(0.7);
  h->GetXaxis()->SetRange(0,500);

  return h;
}


TH1F *getUnbinnedAll(TString &sampleLoc, TString &histoName){
  TH1F *hs = getAndDrawSingleHisto(sampleLoc,"mu0_500","MuonSeed"+histoName);
  hs->SetMarkerSize(0.9);
  hs->SetLineWidth(1);

  hs->GetXaxis()->SetRange(0,500);

  TH1F *hSta = getAndDrawSingleHisto(sampleLoc,"mu0_500","StandAloneMuonsUpdatedAtVtx"+histoName);
  reformatHisto(hSta,22,3);
  hSta->SetMarkerSize(1.3);
  hSta->SetLineWidth(2);
  
  TH1F *hGlb=getAndDrawSingleHisto(sampleLoc,"mu0_500","GlobalMuons"+histoName);
  reformatHisto(hGlb,23,4);
  hGlb->SetMarkerSize(1.3);
  hGlb->SetLineWidth(3);
  
  TH1F *hTk=getAndDrawSingleHisot(sampleLoc,"mu0_500","GeneralTracks"+histoName);
  reformatHisto(hTk,24,5);
  hTk->SetMarkerSize(1.1);
  hTk->SetLineWidth(2);

  return hs;
}
