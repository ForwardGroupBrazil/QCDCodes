
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TAxis.h"
#include "histtools.C"

int plot_two_source_xRay(TString &title="",TString &title2="",TString type){
  gROOT->LoadMacro("adamStyles.C");
  gROOT->LoadMacro("adamCanvasMacros.C");  

  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();
  
  //TString type("PtRate");

  TCanvas *c1 = new TCanvas(title+"/glb_"+type,"Global Muons");
  TCanvas *c2 = new TCanvas(title+"/sta_"+type,"Stand-alone Muons");
  TCanvas *c3 = new TCanvas(title+"/tkmu_"+type,"Tracker Muons");
  
  // Here we insert the cross section and filterEff for the given simulation
  //double x_section_mb = 0.000317; // ttbar xsection
  //double filterEff = 0.33; // ttbar
  
  TFile *f2 = new TFile(title2+".root");
  f2->cd(type);
  TDirectory *thisDir = gDirectory;
  THStack* l3DataNonIsoRate = (THStack*)(thisDir->Get("l3"+type));
  THStack* l2DataNonIsoRate = (THStack*)(thisDir->Get("l2"+type));
  THStack* tkDataNonIsoRate = (THStack*)(thisDir->Get("tkTrack"+type));
  TFile *f1 = new TFile(title+".root");
  //a1 THStack* l3NonIsoRate = (THStack*)thisDir->Get("l3PtRate_cd");
  f1->cd(type);
  thisDir = gDirectory;
  THStack* l3NonIsoRate = (THStack*)thisDir->Get("l3"+type);
  THStack* l2NonIsoRate = (THStack*)thisDir->Get("l2"+type);
  THStack* tkNonIsoRate = (THStack*)thisDir->Get("tkTrack"+type);
  

  c1->cd();
  hist::colors(l3NonIsoRate,1);
  TLegend *leg1 = hist::legend(l3NonIsoRate,"lf");
  leg1->SetHeader("Contribution by #mu parent");
  leg1->SetFillColor(0);
  l3NonIsoRate->Draw("h");
  hist::colors(l3DataNonIsoRate,1);
  l3DataNonIsoRate->Draw("same pe");
  leg1->AddEntry(l3DataNonIsoRate,"Data","p");
  //l3NonIsoRate->GetXaxis()->SetTitle("GLB p_{T} (GeV)");
  l3NonIsoRate->GetYaxis()->SetTitleOffset(1.3);
  float maxY = max(l3NonIsoRate->GetMaximum(),l3DataNonIsoRate->GetMaximum());
  l3NonIsoRate->SetMaximum(1.25*maxY);
  leg1->Draw("same");
  //return 0;

  c2->cd();
  hist::colors(l2NonIsoRate,1);
  TLegend *leg2 = hist::legend(l2NonIsoRate,"lf");
  leg2->SetHeader("Contribution by #mu parent");
  leg2->SetFillColor(0);
  l2NonIsoRate->Draw("h");
  //l2NonIsoRate->GetYaxis()->SetRangeUser(0.,35.);
  hist::colors(l2DataNonIsoRate,1);
  l2DataNonIsoRate->Draw("same pe");
  leg2->AddEntry(l2DataNonIsoRate,"Data","p");
  //l2NonIsoRate->GetXaxis()->SetTitle("STA p_{T} (GeV)");
  l2NonIsoRate->GetYaxis()->SetTitleOffset(1.3);
  maxY = max(l2NonIsoRate->GetMaximum(),l2DataNonIsoRate->GetMaximum());
  l2NonIsoRate->SetMaximum(1.25*maxY);
  leg2->Draw("same");
  
  c3->cd();
  hist::colors(tkNonIsoRate,1);
  TLegend *leg3 = hist::legend(tkNonIsoRate,"lf");
  leg3->SetHeader("Contribution by #mu parent");
  leg3->SetFillColor(0);
  tkNonIsoRate->Draw("h");
  //tkNonIsoRate->GetYaxis()->SetRangeUser(0.,35.);
  hist::colors(tkDataNonIsoRate,1);
  tkDataNonIsoRate->Draw("same pe");
  leg3->AddEntry(tkDataNonIsoRate,"Data","p");
  //tkNonIsoRate->GetXaxis()->SetTitle("TkMu p_{T} (GeV)");
  tkNonIsoRate->GetYaxis()->SetTitleOffset(1.3);
  maxY = max(tkNonIsoRate->GetMaximum(),tkDataNonIsoRate->GetMaximum());
  tkNonIsoRate->SetMaximum(1.25*maxY);
  leg3->Draw("same");

  printCanvasesType(".eps");

  return 0;
  
}
