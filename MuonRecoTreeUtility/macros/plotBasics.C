#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TAxis.h"
#include "histtools.C"

int plotBasics(TString &title="",TString &title2=""){
  gROOT->LoadMacro("adamStyles.C");
  gROOT->LoadMacro("adamCanvasMacros.C");  

  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gStyle->SetErrorX(0.2);
  gStyle->SetEndErrorSize(2);
  gStyle->SetErrorX(0.5);
  gROOT->ForceStyle();

  TFile *f1 = new TFile(title+".root");
  f1->cd("basic_hist");
  TH1* mc_L3_d0 = (TH1*)(gDirectory->Get("hist_L3_d0"));
  TH1* mc_L3_eta = (TH1*)(gDirectory->Get("hist_L3_eta"));
  TH1* mc_L3_pt = (TH1*)(gDirectory->Get("hist_L3_pt"));
  TH1* mc_L3_p = (TH1*)(gDirectory->Get("hist_L3_p"));
  TH1* mc_tk_d0 = (TH1*)(gDirectory->Get("hist_tk_d0"));
  TH1* mc_tk_eta = (TH1*)(gDirectory->Get("hist_tk_eta"));
  TH1* mc_tk_pt = (TH1*)(gDirectory->Get("hist_tk_pt"));
  TH1* mc_tk_p = (TH1*)(gDirectory->Get("hist_tk_p"));


  TFile *f2 = new TFile(title2+".root");
  f2->cd("basic_hist");
  TH1* data_L3_d0 = (TH1*)(gDirectory->Get("hist_L3_d0"));
  TH1* data_L3_eta = (TH1*)(gDirectory->Get("hist_L3_eta"));
  TH1* data_L3_pt = (TH1*)(gDirectory->Get("hist_L3_pt"));
  TH1* data_L3_p = (TH1*)(gDirectory->Get("hist_L3_p"));
  TH1* data_tk_d0 = (TH1*)(gDirectory->Get("hist_tk_d0"));
  TH1* data_tk_eta = (TH1*)(gDirectory->Get("hist_tk_eta"));
  TH1* data_tk_pt = (TH1*)(gDirectory->Get("hist_tk_pt"));
  TH1* data_tk_p = (TH1*)(gDirectory->Get("hist_tk_p"));

  TCanvas *c1 = new TCanvas(title+"/glb_pt","Global Muons");
  c1->cd();
  TH1* mc_L3_pt_clone = mc_L3_pt->Clone("clone");
  mc_L3_pt_clone->SetLineColor(kBlue);
  mc_L3_pt_clone->SetFillColor(kBlue);
  mc_L3_pt_clone->SetFillStyle(3001);
  //mc_L3_pt->Draw("");
  mc_L3_pt_clone->Draw("E2");
  mc_L3_pt_clone->GetYaxis()->SetTitle("muons");
  mc_L3_pt_clone->GetXaxis()->SetTitle("p_{T} [GeV]");
  mc_L3_pt->SetLineColor(kBlue);
  mc_L3_pt->Draw("H same");

  data_L3_pt->SetMarkerStyle(22);
  data_L3_pt->Draw("PE same");

  float maxY = max(mc_L3_pt_clone->GetMaximum(),data_L3_pt->GetMaximum());
  mc_L3_pt_clone->GetYaxis()->SetRangeUser(0.,1.25*maxY);
  
  Float_t xmin = 0.50; 
  Float_t xmax = 0.75;
  Float_t ymin = 0.71;
  Float_t ymax = 0.82;

  TLegend *leg1 = new TLegend(xmin, ymin, xmax, ymax);
  leg1->SetFillColor(0);
  leg1->AddEntry(mc_L3_pt_clone,"MC","lf");
  leg1->AddEntry(data_L3_pt,"Data","lp");
  leg1->Draw("same");

  TCanvas *c2 = new TCanvas(title+"/glb_eta","Global Muons");
  c2->cd();
  TH1* mc_L3_eta_clone = mc_L3_eta->Clone("clone");
  mc_L3_eta_clone->SetLineColor(kBlue);
  mc_L3_eta_clone->SetFillColor(kBlue);
  mc_L3_eta_clone->SetFillStyle(3001);
  //gStyle->SetErrorX(mc_L3_eta->GetBinWidth(3));
  mc_L3_eta_clone->Draw("E2");
  mc_L3_eta_clone->GetYaxis()->SetTitle("muons");
  mc_L3_eta_clone->GetXaxis()->SetTitle("#eta");
  mc_L3_eta->SetLineColor(kBlue);
  mc_L3_eta->Draw("H same");

  maxY = max(mc_L3_eta_clone->GetMaximum(),data_L3_eta->GetMaximum());
  mc_L3_eta_clone->GetYaxis()->SetRangeUser(0.,1.25*maxY);

  data_L3_eta->SetMarkerStyle(22);
  data_L3_eta->Draw("PE same");
  leg1->Draw("same");

  TCanvas *c3 = new TCanvas(title+"/glb_d0","Global Muons");
  c3->cd();
  TH1* mc_L3_d0_clone = mc_L3_d0->Clone("clone");
  mc_L3_d0_clone->SetLineColor(kBlue);
  mc_L3_d0_clone->SetFillColor(kBlue);
  mc_L3_d0_clone->SetFillStyle(3001);
  //gStyle->SetErrorX(mc_L3_d0->GetBinWidth(3));
  mc_L3_d0_clone->Draw("E2");
  mc_L3_d0_clone->GetYaxis()->SetTitle("muons");
  mc_L3_d0_clone->GetXaxis()->SetTitle("d_{0} [cm]");
  mc_L3_d0->SetLineColor(kBlue);
  mc_L3_d0->Draw("H same");

  maxY = max(mc_L3_d0_clone->GetMaximum(),data_L3_d0->GetMaximum());
  mc_L3_d0_clone->GetYaxis()->SetRangeUser(0.,1.25*maxY);

  data_L3_d0->SetMarkerStyle(22);
  data_L3_d0->Draw("PE same");
  //leg1->Draw("same");

  TCanvas *c7 = new TCanvas(title+"/glb_p ","Global Muons");
  c7->cd();
  TH1* mc_L3_p_clone = mc_L3_p->Clone("clone");
  mc_L3_p_clone->SetLineColor(kBlue);
  mc_L3_p_clone->SetFillColor(kBlue);
  mc_L3_p_clone->SetFillStyle(3001);
  //gStyle->SetErrorX(mc_L3_d0->GetBinWidth(3));
  mc_L3_p_clone->Draw("E2");
  mc_L3_p_clone->GetYaxis()->SetTitle("muons");
  mc_L3_p_clone->GetXaxis()->SetTitle("p [GeV]");
  mc_L3_p->SetLineColor(kBlue);
  mc_L3_p->Draw("H same");

  maxY = max(mc_L3_p_clone->GetMaximum(),data_L3_p->GetMaximum());
  mc_L3_p_clone->GetYaxis()->SetRangeUser(0.,1.25*maxY);

  data_L3_p->SetMarkerStyle(22);
  data_L3_p->Draw("PE same");
  //leg1->Draw("same");
  ///////

  TCanvas *c4 = new TCanvas(title+"/tk_pt","Tk Muons");
  c4->cd();
  TH1* mc_tk_pt_clone = mc_tk_pt->Clone("clone");
  mc_tk_pt_clone->SetLineColor(kBlue);
  mc_tk_pt_clone->SetFillColor(kBlue);
  mc_tk_pt_clone->SetFillStyle(3001);
  //gStyle->SetErrorX(mc_tk_pt->GetBinWidth(3));
  mc_tk_pt_clone->Draw("E2");
  mc_tk_pt_clone->GetYaxis()->SetTitle("muons");
  mc_tk_pt_clone->GetXaxis()->SetTitle("p_{T} [GeV]");
  mc_tk_pt->SetLineColor(kBlue);
  mc_tk_pt->Draw("H same");

  maxY = max(mc_tk_pt_clone->GetMaximum(),data_tk_pt->GetMaximum());
  mc_tk_pt_clone->GetYaxis()->SetRangeUser(0.,1.25*maxY);

  data_tk_pt->SetMarkerStyle(22);
  data_tk_pt->Draw("PE same");
  leg1->Draw("same");

  TCanvas *c5 = new TCanvas(title+"/tk_eta","Tk Muons");
  c5->cd();
  TH1* mc_tk_eta_clone = mc_tk_eta->Clone("clone");
  mc_tk_eta_clone->SetLineColor(kBlue);
  mc_tk_eta_clone->SetFillColor(kBlue);
  mc_tk_eta_clone->SetFillStyle(3001);
  //gStyle->SetErrorX(mc_tk_eta->GetBinWidth(3));
  mc_tk_eta_clone->Draw("E2");
  mc_tk_eta_clone->GetYaxis()->SetTitle("muons");
  mc_tk_eta_clone->GetXaxis()->SetTitle("#eta");
  mc_tk_eta->SetLineColor(kBlue);
  mc_tk_eta->Draw("H same");

  maxY = max(mc_tk_eta_clone->GetMaximum(),data_tk_eta->GetMaximum());
  mc_tk_eta_clone->GetYaxis()->SetRangeUser(0.,1.25*maxY);

  data_tk_eta->SetMarkerStyle(22);
  data_tk_eta->Draw("PE same");
  leg1->Draw("same");

  TCanvas *c6 = new TCanvas(title+"/tk_d0","Tk Muons");
  c6->cd();
  TH1* mc_tk_d0_clone = mc_tk_d0->Clone("clone");
  mc_tk_d0_clone->SetLineColor(kBlue);
  mc_tk_d0_clone->SetFillColor(kBlue);
  mc_tk_d0_clone->SetFillStyle(3001);
  //gStyle->SetErrorX(mc_tk_d0->GetBinWidth(3));
  mc_tk_d0_clone->Draw("E2");
  mc_tk_d0_clone->GetYaxis()->SetTitle("muons");
  mc_tk_d0_clone->GetXaxis()->SetTitle("d_{0} [cm]");
  mc_tk_d0->SetLineColor(kBlue);
  mc_tk_d0->Draw("H same");

  maxY = max(mc_tk_d0_clone->GetMaximum(),data_tk_d0->GetMaximum());
  mc_tk_d0_clone->GetYaxis()->SetRangeUser(0.,1.25*maxY);

  data_tk_d0->SetMarkerStyle(22);
  data_tk_d0->Draw("PE same");
  //leg1->Draw("same");

  TCanvas *c8 = new TCanvas(title+"/tk_p ","Tk Muons");
  c8->cd();
  TH1* mc_tk_p_clone = mc_tk_p->Clone("clone");
  mc_tk_p_clone->SetLineColor(kBlue);
  mc_tk_p_clone->SetFillColor(kBlue);
  mc_tk_p_clone->SetFillStyle(3001);
  //gStyle->SetErrorX(mc_tk_d0->GetBinWidth(3));
  mc_tk_p_clone->Draw("E2");
  mc_tk_p_clone->GetYaxis()->SetTitle("muons");
  mc_tk_p_clone->GetXaxis()->SetTitle("p [GeV]");
  mc_tk_p->SetLineColor(kBlue);
  mc_tk_p->Draw("H same");

  maxY = max(mc_tk_p_clone->GetMaximum(),data_tk_p->GetMaximum());
  mc_tk_p_clone->GetYaxis()->SetRangeUser(0.,1.25*maxY);

  data_tk_p->SetMarkerStyle(22);
  data_tk_p->Draw("PE same");
  //leg1->Draw("same");

  printCanvasesType(".png");

  //f1->Close();
  //f2->Close();

  return 0;

}
