#include "TGraph.h"
#include "TGraphErrors.h"
#include <vector>

TH1 *getAndDrawSingleHisto(TString & sampleLoc, TString & fileType, TString & histoName, TString newName = ""){
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

  TGraph * myGraph = drawHisto(histo);

  TString type = reformatLegend(fileType);
  if(newName > "") gLegend->AddEntry(myGraph,newName,"PL");
  else gLegend->AddEntry(myGraph, type, "PL");

  return histo;
}

TH1 *getAndDrawSlicedHistoOrig(TString & sampleLoc, TString & fileType, TString & histoName, int hNum = 2, TString newName = "",float fitRange = 0.0){  
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

  float fitMax = (fitRange>0.) ? fitRange : histo->GetYaxis()->GetXmax();
  float fitMin = (fitRange>0.) ? -1*fitRange : histo->GetYaxis()->GetXmin();

  TF1 *fitFcn = 0;
  fitFcn = new TF1("gauss","gaus",fitMin,fitMax);

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

  TGraph* myGraph = drawHisto(myHisto);

  TString type = reformatLegend(fileType);
  if(newName > "") gLegend->AddEntry(myGraph,newName,"PL");
  else gLegend->AddEntry(myGraph, type, "PL");

  return myHisto;
}


TGraph* drawHisto(TH1 *histo){
  histo->SetLineColor(gLineColor);
  histo->SetMarkerColor(gLineColor);  
  histo->SetMarkerStyle(gMarkerStyle);

  TString optColor(gLineColor==1 ? "" : "same");
  if(gForce) optColor = "same";
  TString optErr = (gOpt==1 ? "E1X0" : "") + optColor;

  histo->Clear();
  histo->Draw(optErr);

  TGraph *myGraph;
  myGraph = ConnectLines(histo);
  if (gLineColor==28) gLineColor=3;
  ++gLineColor;
  ++gLineStyle;
  ++gMarkerStyle;
  if (gLineColor == 5 || gLineColor == 10) ++gLineColor;
  return myGraph;
}


TGraph* ConnectLines(TH1* hist)
{
  TGraph * myGraph = new TGraph(hist);
  myGraph->SetLineColor(gLineColor);
  myGraph->SetLineStyle(gLineStyle);
  myGraph->SetMarkerStyle(hist->GetMarkerStyle());
  if(gOpt==1) myGraph->Draw("L");
  return myGraph;
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

TH1 *getHisto(TString & sampleLoc, TString & fileType, TString & histoName){
  TString fileName = sampleLoc+fileType+".root"; 
  cout<<"Opening the file: "<<fileName<<endl;
  TFile * file =  new TFile(fileName); file->cd();
  if(!file->IsOpen()){
    cout<<"File not opened!"<<endl;
    return;
  }
  cout<<"Getting the file: "<<histoName<<endl;
  TH1 *histo = file->IsOpen() ? (TH1*)file->Get(histoName) : 0;
  if(!histo){ 
    cout<<"Histo not found!"<<endl;  
    return;
  }
  cout<<"Entries: " << histo->GetEntries()<<endl;
  return histo;
}

TH1 *computeEfficiency(const TH1F *num, const TH1F *denom){
  
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

TH1 *getAndDrawSlicedHisto(TString & sampleLoc, TString & fileType, TString & histoName, int hNum = 2, TString newName = "",float fitRange = 0.0){
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

  TGraph* result = fit(histo);      cout << "Line 442" << endl;      
  TGraph * graph= result;
  if(graph)  cout<<"YES!"<<endl;
  drawGraph(graph);
  //graph->Clear();

  TString type = reformatLegend(fileType);
  if(newName > "") gLegend->AddEntry(graph,newName,"PL");
  else gLegend->AddEntry(graph, type, "PL");

  return graph->GetHistogram();
}


TGraph* fit(const TH2* histo, int fitType = 2){
  vector<const TGraphErrors*> results;

  static int i = 0;
  i++;

  vector<double> width;
  vector<double> mass;
  // Errors
  vector<double> widthError;
  vector<double> massError;
  // Chi2
  vector<double> chi2;
  //>>

  //X bin centre and width
  vector<double> binCentre;
  vector<double> binWidth;

  // Fit slices projected along Y from bins in X 
  double cont_min = 100;    //Minimum number of entries
  Int_t binx =  histo->GetXaxis()->GetNbins();

  //  histo->RebinY(4);
  //  histo->RebinX(2);

  for (int i = 1; i < binx ; i++) {
    TH1 *histoY =  histo->ProjectionY(" ", i, i);
    double cont = histoY->GetEntries();
    //cout << "Entries: " << cont << endl;
    if (cont >= cont_min) {
      float minfit = histoY->GetMean() - histoY->GetRMS();
      float maxfit = histoY->GetMean() + histoY->GetRMS();
      TString iString(i);
      TF1 *fitFcn = new TF1(TString("g")+iString,"gaus",minfit,maxfit);
      double x1,x2;
      fitFcn->GetRange(x1,x2);
      //cout << "Range: " << x1 << " " << x2 << endl; 

      histoY->Fit(fitFcn,"0","",x1,x2);

//      histoY->Fit(fitFcn->GetName(),"RME");
      double *par = fitFcn->GetParameters();
      double *err = fitFcn->GetParErrors();

      // Values
      mass.push_back(par[1]);
      width.push_back(par[2]);

      // Errors
      massError.push_back(err[1]);
      widthError.push_back(err[2]);

      chi2.push_back(fitFcn->GetChisquare());
      
      double xx= histo->GetXaxis()->GetBinCenter(i);
      binCentre.push_back(xx);
      
      double ex = 0; //FIXME: you can use the bin width
      binWidth.push_back(ex); 

    }
    else continue;
  }

  // Put the fit results in arrays for TGraphErrors
  const int nn= mass.size();
  double *x = new double[nn];
  double *ym = new double[nn];
  double *e = new double[nn];
  double *eym = new double[nn];
  double *yw = new double[nn];
  double *eyw = new double[nn];
  double *yc = new double[nn];
  
  for (int j=0;j<nn;j++){
    x[j]=binCentre[j];
    ym[j]=mass[j];
    eym[j]=massError[j];
    yw[j]=width[j];
    eyw[j]=widthError[j];
    //     if(j==nn-1)
    //        if(eyw[j] > (eyw[j-1]+0.6)) eyw[j]-=0.6;
    //        if(j==0){
    // 	 if(eyw[0] > (eyw[1]+0.6)) eyw[0]-=0.6;}
    //        else
    // 	 if(eyw[j] > (eyw[j-1]+0.6)) eyw[j]-=0.6;
    //     if(eyw[j] > 1.) eyw[j]=1;
    
    yc[j]=chi2[j];
    e[j]=binWidth[j];
  }

  //Create TGraphErrors
  TString name = histo->GetName();
  TGraphErrors *grM = new TGraphErrors(nn,x,ym,e,eym);
  grM->SetTitle(name+"_Mean");
  grM->SetName(name+"_Mean");
  TGraphErrors *grW = new TGraphErrors(nn,x,yw,e,eyw);
  grW->SetTitle(name+"_Width");
  grW->SetName(name+"_Width");
  TGraphErrors *grC = new TGraphErrors(nn,x,yc,e,e);
  grC->SetTitle(name+"_chi2");
  grC->SetName(name+"_chi2");


  //results.first = grM;  cout << "Line 706" << endl;      
  //results.second = grW;  cout << "Line 707" << endl;      

  //  results.push_back(grC);

  //if(grW)  cout<<"YES! 709"<<endl;

  //  return result; 


  if(fitType==1) return grM;
  if(fitType==2) return grW;
  if(fitType==3) return grC;
  return grW;

}

void drawGraph(TGraph *graph){
  if(!graph) return;
  graph->SetLineColor(gLineColor);
  graph->SetMarkerColor(gLineColor);
  graph->SetMarkerStyle(gMarkerStyle);

  TString optErr(gLineColor==1 ? "PAL" : "PL");

  graph->Clear();
  graph->Draw(optErr);

  if (gLineColor==28) gLineColor=3;
  ++gLineColor;
  ++gLineStyle;
  ++gMarkerStyle;
  if (gLineColor == 5 || gLineColor == 10) ++gLineColor;  

}
