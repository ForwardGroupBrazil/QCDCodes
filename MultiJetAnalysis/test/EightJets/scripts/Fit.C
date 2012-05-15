#include <iostream>
#include <string.h>
#include <fstream>
#include <iomanip>
using namespace RooFit;
void Fit(TString SIGNAL,TString HISTO,double scaleSGN)
{
  gROOT->ForceStyle();
  TFile *fDat = TFile::Open("Histo_flatTree_data_tmva"+SIGNAL+".root");
  TFile *fBkg = TFile::Open("Histo_flatTree_qcd_weights_tmva"+SIGNAL+".root");
  TFile *fSgn = TFile::Open("Histo_flatTree_"+SIGNAL+"_weights_tmva"+SIGNAL+".root");
  
  TH1 *hDat = (TH1*)fDat->Get(HISTO);
  TH1 *hBkgRaw = (TH1*)fBkg->Get(HISTO);
  TH1 *hSgn = (TH1*)fSgn->Get(HISTO);
  TH1 *hDat_JESlo = (TH1*)fDat->Get(HISTO+"_JESlo");
  TH1 *hBkgRaw_JESlo = (TH1*)fBkg->Get(HISTO+"_JESlo");
  TH1 *hSgn_JESlo = (TH1*)fSgn->Get(HISTO+"_JESlo");
  TH1 *hDat_JESup = (TH1*)fDat->Get(HISTO+"_JESup");
  TH1 *hBkgRaw_JESup = (TH1*)fBkg->Get(HISTO+"_JESup");
  TH1 *hSgn_JESup = (TH1*)fSgn->Get(HISTO+"_JESup");
  
  TH1F *hBkg = (TH1F*)hBkgRaw->Clone("Bkg");
  TH1F *hBkg_JESlo = (TH1F*)hBkgRaw_JESlo->Clone("Bkg_JESlo");
  TH1F *hBkg_JESup = (TH1F*)hBkgRaw_JESup->Clone("Bkg_JESup");
  
  hBkg->Smooth(2);
  hBkg_JESlo->Smooth(2);
  hBkg_JESup->Smooth(2);
  hSgn->Smooth(2);
  hSgn_JESlo->Smooth(2);
  hSgn_JESup->Smooth(2);
  
  double lumi = 4967;
  hBkg->Scale(lumi);
  hBkg_JESlo->Scale(lumi);
  hBkg_JESup->Scale(lumi);
  double k_factor = hDat->Integral()/hBkg->Integral();
  double k_factor_JESlo = hDat->Integral()/hBkg_JESlo->Integral();
  double k_factor_JESup = hDat->Integral()/hBkg_JESup->Integral();
  hBkg->Scale(k_factor);
  cout<<"Signal entries = "<<hSgn->GetEntries()<<endl;
  hSgn->Scale(lumi/scaleSGN);
  hBkg_JESlo->Scale(k_factor_JESlo);
  hSgn_JESlo->Scale(lumi/scaleSGN);
  hBkg_JESup->Scale(k_factor_JESup);
  hSgn_JESup->Scale(lumi/scaleSGN);
  hSgn_JESlo->Scale(hSgn->Integral()/hSgn_JESlo->Integral());
  hSgn_JESup->Scale(hSgn->Integral()/hSgn_JESup->Integral());
  
  TH1 *hBkg_STATlo = (TH1*)hBkg->Clone(HISTO+"_STATlo");
  TH1 *hSgn_STATlo = (TH1*)hSgn->Clone(HISTO+"_STATlo");
  TH1 *hBkg_STATup = (TH1*)hBkg->Clone(HISTO+"_STATup");
  TH1 *hSgn_STATup = (TH1*)hSgn->Clone(HISTO+"_STATup");
  
  float y1,e1;
  for(int i=0;i<hBkg->GetNbinsX();i++) {
    y1 = hBkg->GetBinContent(i+1);
    e1 = hBkg->GetBinError(i+1);
    hBkg_STATlo->SetBinContent(i+1,y1-e1);
    hBkg_STATup->SetBinContent(i+1,y1+e1);
    y1 = hSgn->GetBinContent(i+1);
    e1 = hSgn->GetBinError(i+1);
    hSgn_STATlo->SetBinContent(i+1,y1-e1);
    hSgn_STATup->SetBinContent(i+1,y1+e1);
  }
  hBkg_STATlo->Scale(hBkg->Integral()/hBkg_STATlo->Integral());
  hBkg_STATup->Scale(hBkg->Integral()/hBkg_STATup->Integral());
  hSgn_STATlo->Scale(hSgn->Integral()/hSgn_STATlo->Integral());
  hSgn_STATup->Scale(hSgn->Integral()/hSgn_STATup->Integral());
  
  double xMIN = hBkg->GetBinLowEdge(1);
  double xMAX = hBkg->GetBinLowEdge(hBkg->GetNbinsX()+1);
  double xMIN2 = hDat->GetBinLowEdge(hDat->FindFirstBinAbove(0.5));
  double xMAX2 = hDat->GetBinLowEdge(hDat->FindLastBinAbove(0.5)+1);
  RooRealVar x("x","x",xMIN2,xMAX2);
  
  RooDataHist data("data","dataset with x",x,hDat);
  RooDataHist bkg("qcd","bkg with x",x,hBkg);
  RooDataHist sgn("signal","sgn with x",x,hSgn);
  
  RooHistPdf bkgPDF("bkgPDF","bkgPDF",x,bkg,0);
  RooHistPdf sgnPDF("sgnPDF","sgnPDF",x,sgn,0);
  
  RooRealVar f("f","f",0,0.,1.);
  
  RooAddPdf model("model","model",RooArgList(sgnPDF,bkgPDF),RooArgList(f));
  
  RooFitResult* r = model.fitTo(data,Save());
  r->Print("v");
  double N = hDat->Integral();
  double B = hBkg->Integral();
  double S = hSgn->Integral();
  double m = f.getVal();
  double e = f.getError();
  cout<<"k-factor = "<<k_factor<<endl;
  cout<<N<<" "<<B<<" "<<S<<endl;
  cout<<"Total cross section =       "<<N/lumi<<" pb"<<endl;
  cout<<"Model cross section =       "<<S/lumi<<" pb"<<endl;
  cout<<"Fitted signal strength =    "<<m*N/S<<endl;
  cout<<"Fitted signal error =       "<<e*N/S<<endl;
  double p = 0.95;
  double xup = (N/S)*(m+sqrt(2.)*e*TMath::ErfInverse((1-p)*TMath::Erf(m/e)+p));
  cout<<"Bayesian Upper limit =      "<<xup<<endl;
  RooPlot* frame1 = x.frame();  
  data.plotOn(frame1);
  model.plotOn(frame1,Components("sgnPDF*"),LineStyle(1),LineWidth(2),LineColor(kGreen+1));
  model.plotOn(frame1,Components("bkgPDF*"),LineStyle(1),LineWidth(2),LineColor(kRed));
  model.plotOn(frame1,LineStyle(1),LineWidth(2),LineColor(kBlue));
  
  //cout<<frame1->chiSquare()<<endl;
  RooHist* hresid = frame1->residHist();
  RooHist* hpull = frame1->pullHist();
  RooPlot* frame2 = x.frame(Title("Residual Distribution"));
  frame2->addPlotable(hresid,"P") ;
  
  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* frame3 = x.frame(Title("Pull Distribution"));
  frame3->addPlotable(hpull,"P");
  
  TCanvas* cFit = new TCanvas("fitANN_"+SIGNAL,"fitANN_"+SIGNAL,900,600);
  gPad->SetLogy();
  frame1->SetMaximum(1e+4);
  frame1->SetMinimum(0.5);
  frame1->GetXaxis()->SetTitle("ANN Output");
  frame1->GetYaxis()->SetTitle("Events");
  frame1->Draw();
  cout<<frame1->nameOf(3)<<endl;
  TLegend *leg = new TLegend(0.7,0.65,0.9,0.9);
  leg->SetHeader(SIGNAL);
  leg->AddEntry(frame1->findObject("h_data"),"data","P");
  leg->AddEntry(frame1->findObject("model_Norm[x]"),"QCD+Signal","L");
  leg->AddEntry(frame1->findObject("model_Norm[x]_Comp[bkgPDF*]"),"QCD","L");
  leg->AddEntry(frame1->findObject("model_Norm[x]_Comp[sgnPDF*]"),"Signal","L");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->Draw();

  TCanvas* cPull = new TCanvas("pullANN_"+SIGNAL,"pullANN_"+SIGNAL,900,400); 
  frame3->GetXaxis()->SetTitle("ANN Output");
  frame3->GetYaxis()->SetTitle("Pull");
  frame3->Draw();
  
  cout<<"Creating datacard"<<endl;
  ofstream datacard;
  datacard.open("datacard_"+SIGNAL+"_"+HISTO+".txt");
  datacard.setf(ios::right);
  datacard<<"imax 1"<<"\n";
  datacard<<"jmax 1"<<"\n";
  datacard<<"kmax *"<<"\n";
  datacard<<"----------------"<<"\n";
  datacard<<"shapes * * "<<SIGNAL+"_"+HISTO+"_input.root $PROCESS $PROCESS_$SYSTEMATIC"<<"\n";
  datacard<<"----------------"<<"\n";
  datacard<<"bin         1"<<"\n";
  datacard<<"observation "<<N<<"\n";
  datacard<<"----------------"<<"\n";
  datacard<<"bin         1       1"<<"\n";
  datacard<<"process     signal  background  "<<"\n";
  datacard<<"process     0       1"<<"\n";
  datacard<<"rate       "<<S<<"    "<<B<<"\n";
  datacard<<"----------------"<<"\n";
  datacard<<"lumi        lnN       1.022    1.022"<<"\n";
  datacard<<"jes         shape     1        1"<<"\n";
  datacard<<"mcstat      shape     1        1"<<"\n";
  datacard<<"jer         shape     0        0"<<"\n";
  datacard.close();
  
  TFile *out = new TFile(SIGNAL+"_"+HISTO+"_input.root","RECREATE");
  out->cd();
  hDat->Write("data_obs");
  hBkg->Write("background");
  hSgn->Write("signal");
  hDat_JESlo->Write("data_obs_jesDown");
  hBkg_JESlo->Write("background_jesDown");
  hBkg_STATlo->Write("background_mcstatDown");
  hSgn_STATlo->Write("signal_mcstatDown");
  hSgn_JESlo->Write("signal_jesDown");
  hDat_JESup->Write("data_obs_jesUp");
  hBkg_JESup->Write("background_jesUp");
  hBkg_STATup->Write("background_mcstatUp");
  hSgn_JESup->Write("signal_jesUp");
  hSgn_STATup->Write("signal_mcstatUp");
  //----- JER placeholder ----------------
  hBkg_JESlo->Write("background_jerDown");
  hSgn_JESlo->Write("signal_jerDown"); 
  hBkg_JESup->Write("background_jerUp");
  hSgn_JESup->Write("signal_jerUp");
}
