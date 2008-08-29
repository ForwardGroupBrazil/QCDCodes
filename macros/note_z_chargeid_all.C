void note_chargeid_eta(TString &dir="./"){
  gROOT->LoadMacro("adamStyles.C");

  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();

  gROOT->LoadMacro("adamCanvasMacros.C");
  gROOT->LoadMacro("adamDrawMacros.C");
  gROOT->LoadMacro("adamDrawLoop.C");

  gROOT->LoadMacro("chargeId.C");

  gStyle->SetTitleFillColor(0);

  dir = "/scratch/scratch96/a/aeverett/reReco_2011a/postMu/";
  TString dirTkmu =  "/scratch/scratch96/a/aeverett/reReco_2011a/postTkmu/";
  /////
  // MuonSeed _0_
  /////
  TCanvas *c0_01 = newCanvas(dir+"FigParam/resolutions/muonSeed/seed_Idprob_vs_eta","Muon reconstruction efficiences");
  TH1F *h0_01 = chargeId(dir,"MuonSeed","VsEta");
  h0_01->GetXaxis()->SetTitle("#eta");

  TCanvas *c0_02 = newCanvas(dir+"FigParam/resolutions/muonSeed/seed_Idprob_vs_phi","Muon reconstruction efficiences");
  TH1F *h0_02 = chargeId(dir,"MuonSeed","VsPhi");
  h0_02->GetXaxis()->SetTitle("#phi [rad]");

  TCanvas *c1_01 = newCanvas(dir+"FigParam/resolutions/standAlone/sta_Idprob_vs_eta","Muon reconstruction efficiences");
  TH1F *h1_01 = chargeId(dir,"StandAloneMuonsUpdatedAtVtx","VsEta");
  h1_01->GetXaxis()->SetTitle("#eta");

  TCanvas *c1_02 = newCanvas(dir+"FigParam/resolutions/standAlone/sta_Idprob_vs_phi","Muon reconstruction efficiences");
  TH1F *h1_02 = chargeId(dir,"StandAloneMuonsUpdatedAtVtx","VsPhi");
  h1_02->GetXaxis()->SetTitle("#phi [rad]");

  TCanvas *c2_01 = newCanvas(dir+"FigParam/resolutions/tracker/tk_Idprob_vs_eta","Muon reconstruction efficiences");
  TH1F *h2_01 = chargeId(dir,"GeneralTracks","VsEta");
  h2_01->GetXaxis()->SetTitle("#eta");

  TCanvas *c2_02 = newCanvas(dir+"FigParam/resolutions/tracker/tk_Idprob_vs_phi","Muon reconstruction efficiences");
  TH1F *h2_02 = chargeId(dir,"GeneralTracks","VsPhi");
  h2_01->GetXaxis()->SetTitle("#phi [rad]");

  TCanvas *c3_01 = newCanvas(dir+"FigParam/resolutions/global/glb_Idprob_vs_eta","Muon reconstruction efficiences");
  TH1F *h3_01 = chargeId(dir,"GlobalMuons","VsEta");
  h3_01->GetXaxis()->SetTitle("#eta");

  TCanvas *c3_02 = newCanvas(dir+"FigParam/resolutions/global/glb_Idprob_vs_phi","Muon reconstruction efficiences");
  TH1F *h3_02 = chargeId(dir,"GlobalMuons","VsPhi");
  h3_02->GetXaxis()->SetTitle("#phi [rad]");

  TCanvas *c4_01 = newCanvas(dirTkmu+"FigParam/resolutions/trackerMuons/tkmu_Idprob_vs_eta","Muon reconstruction efficiences");
  TH1F *h4_01 = chargeId(dirTkmu,"TrackerMuonsTrackerOnly","VsEta");
  h4_01->GetXaxis()->SetTitle("#eta");

  TCanvas *c4_02 = newCanvas(dirTkmu+"FigParam/resolutions/trackerMuons/tkmu_Idprob_vs_phi","Muon reconstruction efficiences");
  TH1F *h4_02 = chargeId(dirTkmu,"TrackerMuonsTrackerOnly","VsPhi");
  h4_01->GetXaxis()->SetTitle("#phi [rad]");

  //printCanvasesType(".eps");
  //printCanvasesType(".pdf");
  //printCanvasesPS("allPostMuId.ps");
  

}
