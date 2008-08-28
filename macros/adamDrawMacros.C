TH1F *getAndDrawSingleHisto(TString & sampleLoc, TString & fileType, TString & histoName, TString newName = ""){
  TString fileName = sampleLoc+fileType+".root"; 
  if(gLineColor==1) cout<<"Opening the file: "<<fileName<<endl;
  TFile * file =  new TFile(fileName); file->cd();
  if(!file->IsOpen()){
    cout<<"File not opened!"<<endl;
    return;
  }
  if(gLineColor==1) cout<<"Getting the histo: "<<histoName<<endl;
  TH1F *histo = file->IsOpen() ? (TH1F*)file->Get(histoName) : 0;
  if(!histo){ 
    cout<<"Histo not found . . . !"<<endl;  
    return;
  }
  if(gLineColor==1) cout<<"Entries: " << histo->GetEntries()<<endl;

  drawHisto(histo);

  TString type = reformatLegend(fileType);
  if(newName > "") gLegend->AddEntry(histo,newName,"PL");
  else gLegend->AddEntry(histo, type, "PL");

  return histo;
}

TH1 *getAndDrawSlicedHisto(TString & sampleLoc, TString & fileType, TString & histoName, int hNum = 2, TString newName = ""){
  TString fileName = sampleLoc+fileType+".root"; 
  cout<<"Opening the file: "<<fileName<<endl;
  TFile * file =  new TFile(fileName); file->cd();
  if(!file->IsOpen()){
    cout<<"File not opened!"<<endl;
    return;
  }
  cout<<"Getting the histogram: "<<histoName<<endl;
  TH2F *histo = file->IsOpen() ? (TH2F*)file->Get(histoName) : 0;
  if(!histo){ 
    cout<<"Histo not found!"<<endl;  
    return;
  }

  TF1 *fitFcn = 0;
  //  fitFcn = new TF1("gauss","gaus",-0.5,0.5);
  fitFcn = new TF1("gauss","gaus",-0.5,0.5);

  histo->FitSlicesY(fitFcn);
  string name(histo->GetName());
  TH1F* h0 = (TH1*)gDirectory->Get((name+"_0").c_str());
  TH1F* h1 = (TH1*)gDirectory->Get((name+"_1").c_str());
  TH1F* h2 = (TH1*)gDirectory->Get((name+"_2").c_str());
  TH1F* h3 = (TH1*)gDirectory->Get((name+"_chi2").c_str());

  TH1* myHisto = h2;
  
  switch (hNum){
  case 0:{
    myHisto = h0;
    break;
  }
  case 1:{
    myHisto = h1;
    break;
  }
  case 2:{
    myHisto = h2;
    break;
  }
  case 3:{
    myHisto = h3;
    break;
  }
  default: {
    myHisto = h2;
    break;
  }
  }
  
  drawHisto(myHisto);

  TString type = reformatLegend(fileType);
  if(newName > "") gLegend->AddEntry(myHisto,newName,"PL");
  else gLegend->AddEntry(myHisto, type, "PL");

  return myHisto;
}


void drawHisto(TH1 *histo){
  histo->SetLineColor(gLineColor);
  histo->SetMarkerColor(gLineColor);  
  histo->SetMarkerStyle(20);

  TString optColor(gLineColor==1 ? "" : "same");
  if(gForce) optColor = "same";
  TString optErr = (gOpt==1 ? "E1X0" : "") + optColor;
  histo->Clear();

  histo->Draw(optErr);

  if(gOpt==1) ConnectLines(histo);
  if (gLineColor==28) gLineColor=3;
  ++gLineColor;
  ++gLineStyle;
  if (gLineColor == 5 || gLineColor == 10) ++gLineColor;
}



void ConnectLines(TH1* hist)
{
  if ( hist->GetEntries() == 0 ) return;
  const int nBin = hist->GetNbinsX();
  TGraph* grp = new TGraph(nBin);

  for(int i=0; i<nBin; i++) {
    double val = hist->GetBinContent(i+1);
    //    if(val<=0) val = 0.0001; 
    grp->SetPoint(i, hist->GetBinCenter(i+1), val);
  }
  grp->SetLineWidth(1);
  grp->SetLineColor(gLineColor);
  grp->SetLineStyle(gLineStyle);//++
  grp->Draw("LSame");
}


TString reformatLegend(TString &fileType){
  
  TString type = "#mu";

  if(fileType=="muminus5") 
    type=type+"^{-} p_{T} = 5 GeV/c";
  else if(fileType=="muplus5")
    type=type+"^{+} p_{T} = 5 GeV/c";
  else if(fileType=="mu5")
    type=type+" p_{T} = 5 GeV/c";

  else if(fileType=="muminus10") 
    type=type+"^{-} p_{T} = 10 GeV/c";
  else if(fileType=="muplus10")
    type=type+"^{+} p_{T} = 10 GeV/c";
  else if(fileType=="mu10")
    type=type+" p_{T} = 10 GeV/c";

  else if(fileType=="muminus50") 
    type=type+"^{-} p_{T} = 50 GeV/c";
  else if(fileType=="muplus50")
    type=type+"^{+} p_{T} = 50 GeV/c";
  else if(fileType=="mu50")
    type=type+" p_{T} = 50 GeV/c";

  else if(fileType=="muminus100") 
    type=type+"^{-} p_{T} = 100 GeV/c";
  else if(fileType=="muplus100")
    type=type+"^{+} p_{T} = 100 GeV/c";
  else if(fileType=="mu100")
    type=type+" p_{T} = 100 GeV/c";

  else if(fileType=="muminus200") 
    type=type+"^{-} p_{T} = 200 GeV/c";
  else if(fileType=="muplus200")
    type=type+"^{+} p_{T} = 200 GeV/c";
  else if(fileType=="mu200")
    type=type+" p_{T} = 200 GeV/c";

  else if(fileType=="muminus500") 
    type=type+"^{-} p_{T} = 500 GeV/c";
  else if(fileType=="muplus500")
    type=type+"^{+} p_{T} = 500 GeV/c";
  else if(fileType=="mu500")
    type=type+" p_{T} = 500 GeV/c";

  else if(fileType=="muminus1000") 
    type=type+"^{-} p_{T} = 1000 GeV/c";
  else if(fileType=="muplus1000")
    type=type+"^{+} p_{T} = 1000 GeV/c";
  else if(fileType=="mu1000")
    type=type+" p_{T} = 1000 GeV/c";

  else if(fileType=="muminus2000") 
    type=type+"^{-} p_{T} = 2000 GeV/c";
  else if(fileType=="muplus2000")
    type=type+"^{+} p_{T} = 2000 GeV/c";
  else if(fileType=="mu2000")
    type=type+" p_{T} = 2000 GeV/c";

  else if(fileType=="muminus3000") 
    type=type+"^{-} p_{T} = 3000 GeV/c";
  else if(fileType=="muplus3000")
    type=type+"^{+} p_{T} = 3000 GeV/c";
  else if(fileType=="mu3000")
    type=type+" p_{T} = 3000 GeV/c";

  else if(fileType=="muminus1_200") 
    type=type+"^{-} p_{T} = 1-200 GeV/c";
  else if(fileType=="muplus1_200")
    type=type+"^{+} p_{T} = 1-200 GeV/c";
  else if(fileType=="mu1_200")
    type=type+" p_{T} = 1-200 GeV/c";

  else if(fileType=="muminus0_500") 
    type=type+"^{-} p_{T} = 0-500 GeV/c";
  else if(fileType=="muplus0_500")
    type=type+"^{+} p_{T} = 0-500 GeV/c";
  else if(fileType=="mu0_500")
    type=type+" p_{T} = 0-500 GeV/c";

  else if(fileType=="StandAloneMuonsUpdatedAtVtx")
    type="Stand Alone Muons Z";
  else if(fileType=="Muons")
    type="Global Muons Z";

  else
    cout<<"### Case not handled! ###"<<endl;
  
  return type;
}

void reformatHisto(TH1 *histo, int msty, int lsty){
  histo->SetLineStyle(gLineStyle-1);//lsty
  histo->SetMarkerStyle(msty);
  if(msty!=21) histo->SetMarkerSize(1.1);
  TH1F *h2 = histo->Clone();
  h2->SetLineStyle(1);
  h2->Draw("E1X0 same");
}
