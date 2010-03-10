
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TAxis.h"
#include "histtools.C"

int plot_xRay(TString &title=""){
  gROOT->LoadMacro("adamStyles.C");
  
  setHLTStyle();
  gROOT->SetStyle("hltStyle");
  gROOT->ForceStyle();
  
  
  TCanvas *c1 = new TCanvas("glb","Global Muons");
  TCanvas *c2 = new TCanvas("sta","Stand-alone Muons");
  TCanvas *c3 = new TCanvas("tkmu","Tracker Muons");
  
  // Here we insert the cross section and filterEff for the given simulation
  //double x_section_mb = 0.000317; // ttbar xsection
  //double filterEff = 0.33; // ttbar
  
  TFile f1(title+".root");
  //a1 THStack* l3NonIsoRate = (THStack*)rate_hist.Get("l3PtRate_cd");
  THStack* l3NonIsoRate = (THStack*)rate_hist.Get("l3PtRate");
  THStack* l2NonIsoRate = (THStack*)rate_hist.Get("l2PtRate");
  THStack* tkNonIsoRate = (THStack*)rate_hist.Get("tkTrackPtRate");
  

  c1->cd();
  hist::colors(l3NonIsoRate,1);
  TLegend *leg1 = hist::legend(l3NonIsoRate);
  leg1->SetHeader("Contribution by #mu parent");
  l3NonIsoRate->Draw("p");
  l3NonIsoRate->GetXaxis()->SetTitle("GLB p_{T} (GeV)");

  l3NonIsoRate->GetYaxis()->SetTitleOffset(1.3);
  leg1->Draw("same");
  //return 0;
  c2->cd();
  hist::colors(l2NonIsoRate,1);
  TLegend *leg2 = hist::legend(l2NonIsoRate);
  leg2->SetHeader("Contribution by #mu parent");
  l2NonIsoRate->Draw("p");
  l2NonIsoRate->GetXaxis()->SetTitle("STA p_{T} (GeV)");

  l2NonIsoRate->GetYaxis()->SetTitleOffset(1.3);
  leg2->Draw("same");
  
  c3->cd();
  hist::colors(tkNonIsoRate,1);
  TLegend *leg3 = hist::legend(tkNonIsoRate);
  leg3->SetHeader("Contribution by #mu parent");
  tkNonIsoRate->Draw("p");
  tkNonIsoRate->GetXaxis()->SetTitle("TkMu p_{T} (GeV)");

  tkNonIsoRate->GetYaxis()->SetTitleOffset(1.3);
  leg3->Draw("same");

  return 0;
  
}
