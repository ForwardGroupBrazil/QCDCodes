void muonZEfficiencies(TString &dir="./"){
  gROOT->LoadMacro("adamStyles.C");
  gROOT->LoadMacro("adamCanvasMacros.C");
  gROOT->LoadMacro("adamDrawMacros.C");
  gROOT->LoadMacro("adamDrawLoop.C");

  setRBellanStyle();
  gROOT->SetStyle("rbStyle");
  gROOT->ForceStyle();

  gStyle->SetTitleFillColor(0);

  dir = "./";

  /////
  // MuonSeed _0_
  /////
  TCanvas *c0_01 = newCanvas("seed_sim_eff_vs_eta","Muon reconstruction efficiences");
  TH1F *h0_01 = adamDrawLoop(dir,"MuonSeedSimEffVsEta",0,1);

  TCanvas *c0_02 = newCanvas("seed_sim_eff_vs_phi","Muon reconstruction efficiences");
  TH1F *h0_02 = adamDrawLoop(dir,"MuonSeedSimEffVsPhi",0,1);
  h0_02->GetXaxis()->SetTitle("#phi (rad)");
  ///////////////
  ///////////////

  gStyle->SetOptLogy(1);

  TCanvas *c0_03 = newCanvas("seed_res_vs_eta","Muon reconstruction efficiences",0,0,2);
  TH1F *h0_03 = adamDrawLoop(dir,"MuonSeedResoQOverPtVsEta",0,1,2);
  h0_03->GetYaxis()->SetTitle("Resolution");

  TCanvas *c0_04 = newCanvas("seed_res_vs_phi","Muon reconstruction efficiences",0,0,2);
  TH1F *h0_04 = adamDrawLoop(dir,"MuonSeedResoQOverPtVsPhi",0,1,2);
  h0_04->GetYaxis()->SetTitle("Resolution");
  h0_04->GetXaxis()->SetTitle("#phi (rad)");

  printCanvasesType(".eps");
  printCanvasesPS("all.ps");
  
}
