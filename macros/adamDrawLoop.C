//
TLegend *gLegend;

int gLineColor;
int gLineStyle;
int gMarkerStyle;

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
  gMarkerStyle=21;

  gOpt = opt;
  gType = type;

  
  switch(type){
  case 1: {
    //gLegend = new TLegend(.17, .69, 0.39, .92, "");
    gLegend = new TLegend(.35, .13, 0.65, .33, "");
    break;
  }
  case 2: {
    gLegend = new TLegend(.35, .13, 0.65, .33, "");
    //gLegend = new TLegend(.15, .76, 0.32, .93, "");
    break;
  }
  default: {
    gLegend = new TLegend(.67, .10, 0.84, .27, "");
    break;
  }
  }
  

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
    histo = (gType==1) ? getAndDrawSingleHisto(sampleLoc,"mu5",histoName) :  getAndDrawSlicedHisto(sampleLoc,"mu5",histoName);
    break;
  }
  case 10:{
    histo = (gType==1) ? getAndDrawSingleHisto(sampleLoc,"mu10",histoName) : getAndDrawSlicedHisto(sampleLoc,"mu10",histoName);
    break;
  }
  case 50:{
    histo = (gType==1) ? getAndDrawSingleHisto(sampleLoc,"mu50",histoName) : getAndDrawSlicedHisto(sampleLoc,"mu50",histoName);
    break;
  }
  case 100:{
    if(gType==1) histo = getAndDrawSingleHisto(sampleLoc,"mu100",histoName);
    if(gType==2) histo = getAndDrawSlicedHisto(sampleLoc,"mu100",histoName);
    if(gType==3) histo = getAndDrawSlicedHisto(sampleLoc,"mu100",histoName);
    break;
  }
  case 200:{
    histo = (gType==1) ? getAndDrawSingleHisto(sampleLoc,"mu200",histoName) : getAndDrawSlicedHisto(sampleLoc,"mu200",histoName);
    break;
  }
  case 500:{
    histo = (gType==1) ? getAndDrawSingleHisto(sampleLoc,"mu500",histoName) :  getAndDrawSlicedHisto(sampleLoc,"mu500",histoName);
    break;
  }
  case 1000:{
    if(gType==1) histo = getAndDrawSingleHisto(sampleLoc,"mu1000",histoName);
    if(gType==2) histo = getAndDrawSlicedHisto(sampleLoc,"mu1000",histoName);
    if(gType==3) histo = getAndDrawSlicedHisto(sampleLoc,"mu1000",histoName);
    break;
  }
  case 2000:{
    histo = (gType==1) ? getAndDrawSingleHisto(sampleLoc,"mu2000",histoName) : getAndDrawSlicedHisto(sampleLoc,"mu2000",histoName);
    break;
  }
  case 3000:{
    histo = (gType==1) ? getAndDrawSingleHisto(sampleLoc,"mu3000",histoName) : getAndDrawSlicedHisto(sampleLoc,"mu3000",histoName);
    break;
  }
  case 999:{
    if(gType==1) histo = getAndDrawSingleHisto(sampleLoc,"mu0_500",histoName);
    if(gType==2) histo = getAndDrawSlicedHisto(sampleLoc,"mu0_500",histoName,2);
    if(gType==3) histo = getAndDrawSlicedHisto(sampleLoc,"mu0_500",histoName,1);
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

  switch(xAxis){
  case 0:{
    histo->GetXaxis()->SetTitle("");
    break;
  }
  case 1:{
    histo->GetXaxis()->SetTitle("#eta");
    break;
  }
  case 2:{
    histo->GetXaxis()->SetTitle("#phi (rad)");
    break;
  }
  case 3:{
    histo->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    break;
  }
  }

  switch(type){
  case 1:{
    histo->GetYaxis()->SetTitle("Efficiency");
    histo->SetMinimum(0.9);
    histo->SetMaximum(1.002);
    break;
  }
  case 2:{
    histo->GetYaxis()->SetTitle("Resolution");
    //histo->SetMinimum(0.001);
    //histo->SetMaximum(10.);
    break;
  }
  case 3:{
    histo->GetYaxis()->SetTitle("Pull");
    //histo->SetMinimum(0.);
    ////histo->SetMaximum(2.);
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
  //  if(gType==3) h5 = getAndDrawSlicedHisto(sampleLoc,"mu5",histoName,1);

  TH1 *h10;
  if(gType==1) h10 = getAndDrawSingleHisto(sampleLoc,"mu10",histoName);
  if(gType==2) h10 = getAndDrawSlicedHisto(sampleLoc,"mu10",histoName,2);
  if(gType==3) h10 = getAndDrawSlicedHisto(sampleLoc,"mu10",histoName,1);

  TH1 *h100;
  if(gType==1) h100 = getAndDrawSingleHisto(sampleLoc,"mu100",histoName);
  if(gType==2) h100 = getAndDrawSlicedHisto(sampleLoc,"mu100",histoName,2);
  if(gType==3) h100 = getAndDrawSlicedHisto(sampleLoc,"mu100",histoName,1);

  TH1 *h200;
  if(gType==1) h200 = getAndDrawSingleHisto(sampleLoc,"mu200",histoName);
  if(gType==2) h200 = getAndDrawSlicedHisto(sampleLoc,"mu200",histoName,2);
  if(gType==3) h200 = getAndDrawSlicedHisto(sampleLoc,"mu200",histoName,1);

  TH1 *h500;
  if(gType==1) h500 = getAndDrawSingleHisto(sampleLoc,"mu500",histoName);
  if(gType==2) h500 = getAndDrawSlicedHisto(sampleLoc,"mu500",histoName,2);
  if(gType==3) h500 = getAndDrawSlicedHisto(sampleLoc,"mu500",histoName,1);

  TH1 *h1000;
  if(gType==1) h1000 = getAndDrawSingleHisto(sampleLoc,"mu1000",histoName);
  if(gType==2) h1000 = getAndDrawSlicedHisto(sampleLoc,"mu1000",histoName,2);
  if(gType==3) h1000 = getAndDrawSlicedHisto(sampleLoc,"mu1000",histoName,1);

  /*  
      TH1 *h2000;
      if(gType==1) h2000 = getAndDrawSingleHisto(sampleLoc,"mu2000",histoName);
      if(gType==2) h2000 = getAndDrawSlicedHisto(sampleLoc,"mu2000",histoName,2);
      if(gType==3) h2000 = getAndDrawSlicedHisto(sampleLoc,"mu2000",histoName,1);
      
      TH1 *h3000;
      if(gType==1) h3000 = getAndDrawSingleHisto(sampleLoc,"mu3000",histoName);
      if(gType==2) h3000 = getAndDrawSlicedHisto(sampleLoc,"mu3000",histoName,2);
      if(gType==3) h3000 = getAndDrawSlicedHisto(sampleLoc,"mu3000",histoName,1);
  */
  
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
  if(gType==3) hTk = getAndDrawSlicedHisto(sampleLoc,"mu0_500","GeneralTracks"+histoName,1,"General Tracks");
  hTk->GetXaxis()->SetRange(0,500);

  TH1 *hs;
  if(gType==1) hs = getAndDrawSingleHisto(sampleLoc,"mu0_500","MuonSeed"+histoName,"Muon Seed");
  if(gType==2) hs = getAndDrawSlicedHisto(sampleLoc,"mu0_500","MuonSeed"+histoName,2,"Muon Seed");
  if(gType==3) hs = getAndDrawSlicedHisto(sampleLoc,"mu0_500","MuonSeed"+histoName,1,"Muon Seed");

  TH1 *hSta;
  if(gType==1) hSta = getAndDrawSingleHisto(sampleLoc,"mu0_500","StandAloneMuonsUpdatedAtVtx"+histoName,"Stand Alone Muons");
  if(gType==2) hSta = getAndDrawSlicedHisto(sampleLoc,"mu0_500","StandAloneMuonsUpdatedAtVtx"+histoName,2,"Stand Alone Muons");
  if(gType==3) hSta = getAndDrawSlicedHisto(sampleLoc,"mu0_500","StandAloneMuonsUpdatedAtVtx"+histoName,1,"Stand Alone Muons");
  
  TH1 *hGlb;
  if(gType==1) hGlb = getAndDrawSingleHisto(sampleLoc,"mu0_500","GlobalMuons"+histoName,"Global Muons");
  if(gType==2) hGlb = getAndDrawSlicedHisto(sampleLoc,"mu0_500","GlobalMuons"+histoName,2,"Global Muons");
  if(gType==3) hGlb = getAndDrawSlicedHisto(sampleLoc,"mu0_500","GlobalMuons"+histoName,1,"Global Muons");

  TString tkmuDir = "/scratch/scratch96/a/aeverett/reReco_2011a/postTkmu/";
  //GeneralTracks should be TrackerMuonsTrackerOnly
  TH1 *hTkmu;
  if(gType==1) hTkmu = getAndDrawSingleHisto(tkmuDir,"mu0_500","GeneralTracks"+histoName,"Tracker Muons");
  if(gType==2) hTkmu = getAndDrawSlicedHisto(tkmuDir,"mu0_500","GeneralTracks"+histoName,2,"Tracker Muons");
  if(gType==3) hTkmu = getAndDrawSlicedHisto(tkmuDir,"mu0_500","GeneralTracks"+histoName,1,"Tracker Muons");

  return hTk;
}
