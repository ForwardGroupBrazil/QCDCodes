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
 *
 * xAxis = 1: eta
 * xAxis = 2: phi
 * xAxis = 3: pt
 */

TH1 *adamDrawLoop(TString & sampleLoc,TString & histoName, int iSwitch =0, int opt = 0, int type = 1, int xAxis = 1){
  gLineColor=1;
  gLineStyle=1;

  gOpt = opt;
  gType = type;

  switch(type){
  case 1: {
    gLegend = new TLegend(.57, .15, 0.79, .38, ""); // eff 5 GeV seed
    break;
  }
  case 2: {
    gLegend = new TLegend(.4, .7, 0.58, .88, ""); // eff 5 GeV seed
    break;
  }
  default: {
    gLegend = new TLegend(.67, .10, 0.84, .27, ""); // eff 5 GeV seed
    break;
  }
  }
  //  gLegend = new TLegend(.60, .25, 0.83, .45, ""); // eff 5 GeV seed
  // gLegend = new TLegend(.65, .25, 0.82, .42, "");
  //  gLegend = new TLegend(.60, .25, 0.90, .60, ""); // eff 5 GeV seed

  gLegend->SetFillColor(0);
  
  TH1 *histo = 0;

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
    if(gType==1) histo = getAndDrawSingleHisto(sampleLoc,"mu100",histoName);
    if(gType==2) histo = getAndDrawSlicedHisto(sampleLoc,"mu100",histoName);
    if(gType==3) histo = getAndDrawSlicedHisto(sampleLoc,"mu100",histoName);
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
    if(gType==1) histo = getAndDrawSingleHisto(sampleLoc,"mu1000",histoName);
    if(gType==2) histo = getAndDrawSlicedHisto(sampleLoc,"mu1000",histoName);
    if(gType==3) histo = getAndDrawSlicedHisto(sampleLoc,"mu1000",histoName);
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
    if(gType==1) histo = getAndDrawSingleHisto(sampleLoc,"mu0_500",histoName);
    if(gType==2) histo = getAndDrawSlicedHisto(sampleLoc,"mu0_500",histoName,2);
    if(gType==3) histo = getAndDrawSlicedHisto(sampleLoc,"mu0_500",histoName,0);
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

  //aaa  histo->GetXaxis()->SetTitleSize(0.055); // eta
  //aaa  histo->GetXaxis()->SetTitleOffset(0.55); //eta

  //aaa  histo->GetYaxis()->SetTitleOffset(1.25);
  //aaa  histo->GetYaxis()->SetNdivisions(511);
  //aaa  histo->GetYaxis()->SetLabelSize(0.04);

  switch(xAxis){
  case 1:{
    histo->GetXaxis()->SetTitle("#eta");
    break;
  }
  case 2:{
    histo->GetXaxis()->SetTitle("#phi (rad)");
    //aaa    histo->GetXaxis()->SetTitleSize(0.05); // phi
    //aaa    histo->GetXaxis()->SetTitleOffset(0.75); //phi
    break;
  }
  case 3:{
    histo->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    //aaa    histo->GetXaxis()->SetTitleSize(0.05); // phi
    //aaa    histo->GetXaxis()->SetTitleOffset(0.75); //phi
    break;
  }
  }

  switch(type){
  case 1:{
    histo->GetYaxis()->SetTitle("Efficiency");
    histo->SetMinimum(0.9);
    histo->SetMaximum(1.005);
    break;
  }
  case 2:{
    histo->GetYaxis()->SetTitle("Resolution");
    histo->SetMinimum(0.001);
    histo->SetMaximum(10.);
    break;
  }
  case 3:{
    histo->GetYaxis()->SetTitle("Pull");
    histo->SetMinimum(0.);
    //histo->SetMaximum(2.);
    break;
  }
  }

  histo->Clear();

  gLegend->Draw("same");
  return histo;
}

TH1 *getInclusive(TString &sampleLoc, TString &histoName){
  //  TH1F *h5 = getAndDrawSingleHisto(sampleLoc,"mu5",histoName);
  //  h5->SetMarkerSize(0.7);

  //  TH1 *h5;
  //  if(gType==1) h5 = getAndDrawSingleHisto(sampleLoc,"mu5",histoName);
  //  if(gType==2) h5 = getAndDrawSlicedHisto(sampleLoc,"mu5",histoName,2);
  //  if(gType==3) h5 = getAndDrawSlicedHisto(sampleLoc,"mu5",histoName,0);
  //reformatHisto(h5,21,2);


  TH1 *h10;
  if(gType==1) h10 = getAndDrawSingleHisto(sampleLoc,"mu10",histoName);
  if(gType==2) h10 = getAndDrawSlicedHisto(sampleLoc,"mu10",histoName,2);
  if(gType==3) h10 = getAndDrawSlicedHisto(sampleLoc,"mu10",histoName,2);
  reformatHisto(h10,21,2);

  TH1 *h100;
  if(gType==1) h100 = getAndDrawSingleHisto(sampleLoc,"mu100",histoName);
  if(gType==2) h100 = getAndDrawSlicedHisto(sampleLoc,"mu100",histoName,2);
  if(gType==3) h100 = getAndDrawSlicedHisto(sampleLoc,"mu100",histoName,2);
  reformatHisto(h100,22,3);

  TH1 *h200;
  if(gType==1) h200 = getAndDrawSingleHisto(sampleLoc,"mu200",histoName);
  if(gType==2) h200 = getAndDrawSlicedHisto(sampleLoc,"mu200",histoName,2);
  if(gType==3) h200 = getAndDrawSlicedHisto(sampleLoc,"mu200",histoName,2);
  reformatHisto(h200,23,4);

  TH1 *h10;
  if(gType==1) h500 = getAndDrawSingleHisto(sampleLoc,"mu500",histoName);
  if(gType==2) h500 = getAndDrawSlicedHisto(sampleLoc,"mu500",histoName,2);
  if(gType==3) h500 = getAndDrawSlicedHisto(sampleLoc,"mu500",histoName,2);
  reformatHisto(h500,24,5);

  TH1 *h1000;
  if(gType==1) h1000 = getAndDrawSingleHisto(sampleLoc,"mu1000",histoName);
  if(gType==2) h1000 = getAndDrawSlicedHisto(sampleLoc,"mu1000",histoName,2);
  if(gType==3) h1000 = getAndDrawSlicedHisto(sampleLoc,"mu1000",histoName,2);
  reformatHisto(h1000,26,6);

    TH1 *h2000;
    if(gType==1) h2000 = getAndDrawSingleHisto(sampleLoc,"mu2000",histoName);
    if(gType==2) h2000 = getAndDrawSlicedHisto(sampleLoc,"mu2000",histoName,2);
    if(gType==3) h2000 = getAndDrawSlicedHisto(sampleLoc,"mu2000",histoName,0);
    reformatHisto(h2000,26,7);

    TH1 *h3000;
    if(gType==1) h3000 = getAndDrawSingleHisto(sampleLoc,"mu3000",histoName);
    if(gType==2) h3000 = getAndDrawSlicedHisto(sampleLoc,"mu3000",histoName,2);
    if(gType==3) h3000 = getAndDrawSlicedHisto(sampleLoc,"mu3000",histoName,0);
    reformatHisto(h3000,27,8);


  return h10;
}


TH1 *getUnbinned(TString &sampleLoc, TString &histoName){
  TH1 *h = getAndDrawSingleHisto(sampleLoc,"mu0_500",histoName);
  h->SetMarkerSize(0.7);
  h->GetXaxis()->SetRange(0,500);

  return h;
}


TH1 *getUnbinnedAll(TString &sampleLoc, TString &histoName){

  TH1 *hTk;
  if(gType==1) hTk = getAndDrawSingleHisto(sampleLoc,"mu0_500","GeneralTracks"+histoName,"General Tracks");
  if(gType==2) hTk = getAndDrawSlicedHisto(sampleLoc,"mu0_500","GeneralTracks"+histoName,2,"General Tracks");
  if(gType==3) hTk = getAndDrawSlicedHisto(sampleLoc,"mu0_500","GeneralTracks"+histoName,0,"General Tracks");
  reformatHisto(hTk,21,2);
  //hTk->SetMarkerSize(1.1);
  //hTk->SetLineWidth(2);
  hTk->GetXaxis()->SetRange(0,500);

  TH1 *hs;
  if(gType==1) hs = getAndDrawSingleHisto(sampleLoc,"mu0_500","MuonSeed"+histoName,"Muon Seed");
  if(gType==2) hs = getAndDrawSlicedHisto(sampleLoc,"mu0_500","MuonSeed"+histoName,2,"Muon Seed");
  if(gType==3) hs = getAndDrawSlicedHisto(sampleLoc,"mu0_500","MuonSeed"+histoName,0,"Muon Seed");
  reformatHisto(hs,22,3);
  //hs->SetMarkerSize(0.9);
  //hs->SetLineWidth(1);

  TH1 *hSta;
  if(gType==1) hSta = getAndDrawSingleHisto(sampleLoc,"mu0_500","StandAloneMuonsUpdatedAtVtx"+histoName,"Stand Alone Muons");
  if(gType==2) hSta = getAndDrawSlicedHisto(sampleLoc,"mu0_500","StandAloneMuonsUpdatedAtVtx"+histoName,2,"Stand Alone Muons");
  if(gType==3) hSta = getAndDrawSlicedHisto(sampleLoc,"mu0_500","StandAloneMuonsUpdatedAtVtx"+histoName,0,"Stand Alone Muons");
  reformatHisto(hSta,22,3);
  //hSta->SetMarkerSize(1.3);
  //hSta->SetLineWidth(2);
  
  TH1 *hGlb;
  if(gType==1) hGlb = getAndDrawSingleHisto(sampleLoc,"mu0_500","GlobalMuons"+histoName,"Global Muons");
  if(gType==2) hGlb = getAndDrawSlicedHisto(sampleLoc,"mu0_500","GlobalMuons"+histoName,2,"Global Muons");
  if(gType==3) hGlb = getAndDrawSlicedHisto(sampleLoc,"mu0_500","GlobalMuons"+histoName,0,"Global Muons");
  reformatHisto(hGlb,23,4);
  //hGlb->SetMarkerSize(1.3);
  //hGlb->SetLineWidth(3);

  TString tkmuDir = "/scratch/scratch96/a/aeverett/reReco_2011a/postTkmu/";
  //GeneralTracks should be TrackerMuonsTrackerOnly

  TH1 *hTkmu;
  if(gType==1) hTkmu = getAndDrawSingleHisto(tkmuDir,"mu0_500","GeneralTracks"+histoName,"Tracker Muons");
  if(gType==2) hTkmu = getAndDrawSlicedHisto(tkmuDir,"mu0_500","GeneralTracks"+histoName,2,"Tracker Muons");
  if(gType==3) hTkmu = getAndDrawSlicedHisto(tkmuDir,"mu0_500","GeneralTracks"+histoName,0,"Tracker Muons");
  reformatHisto(hTkmu,25,6);
  //hTkmu->SetMarkerSize(1.1);
  //hTkmu->SetLineWidth(2);

  return hTk;
}
