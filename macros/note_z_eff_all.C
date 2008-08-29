void note_z_eff_all(TString &dir="./"){
  gROOT->LoadMacro("adamStyles.C");

  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();

  gROOT->LoadMacro("adamCanvasMacros.C");
  gROOT->LoadMacro("adamDrawMacros.C");
  gROOT->LoadMacro("adamDrawLoop.C");

  gStyle->SetTitleFillColor(0);

  //TLegend should be     gLegend = new TLegend(.35, .13, 0.65, .33, "");

  dir = "/scratch/scratch96/a/aeverett/reReco_2011a/postMu/";
  TString dirTkmu = "/scratch/scratch96/a/aeverett/reReco_2011a/postTkmu/";

  /////
  // MuonSeed _0_
  /////
  TCanvas *c0_01 = newCanvas(dir+"FigSTA/efficiencies/muonSeed/seed_sim_eff_vs_eta","Muon reconstruction efficiences");
  TH1F *h0_01 = adamDrawLoop(dir,"MuonSeedSimEffVsEta",0,1,1,1);

  TCanvas *c0_02 = newCanvas(dir+"FigSTA/efficiencies/muonSeed/seed_sim_eff_vs_phi","Muon reconstruction efficiences");
  TH1F *h0_02 = adamDrawLoop(dir,"MuonSeedSimEffVsPhi",0,1,1,2);

  //  TCanvas *c0_02 = newCanvas("FigSTA/efficiencies/muonSeed/seed_eff_vs_simhits","Muon reconstruction efficiences");
  //  TH1F *h0_02 = adamDrawLoop(dir,"MuonSeedSimEffVsSimHits",0,1);

  /////
  // StandAloneMuon _1_
  /////
  TCanvas *c1_01 = newCanvas(dir+"FigSTA/efficiencies/standAlone/sta_sim_eff_vs_eta","Muon reconstruction efficiences");
  TH1F *h1_01 = adamDrawLoop(dir,"StandAloneMuonsUpdatedAtVtxSimEffVsEta",0,1,1,1);

  TCanvas *c1_02 = newCanvas(dir+"FigSTA/efficiencies/standAlone/sta_sim_eff_vs_phi","Muon reconstruction efficiences");  
  TH1F *h1_02 = adamDrawLoop(dir,"StandAloneMuonsUpdatedAtVtxSimEffVsPhi",0,1,1,2);

  TCanvas *c1_03 = newCanvas(dir+"FigSTA/efficiencies/standAlone/sta_seed_eff_vs_eta","Muon reconstruction efficiences");  
  TH1F *h1_03 = adamDrawLoop(dir,"StandAloneMuonsUpdatedAtVtxMuonSeedEffVsEta",0,1,1,1);

  TCanvas *c1_04 = newCanvas(dir+"FigSTA/efficiencies/standAlone/sta_seed_eff_vs_phi","Muon reconstruction efficiences");  
  TH1F *h1_04 = adamDrawLoop(dir,"StandAloneMuonsUpdatedAtVtxMuonSeedEffVsPhi",0,1,1,2);

  //  TCanvas *c1_05 = newCanvas("FigSTA/efficiencies/standAlone/sta_seed_eff_vs_simhits","Muon reconstruction efficiences");  
  //  TH1F *h1_05 = adamDrawLoop(dir,"StandAloneMuonsUpdatedAtVtxMuonSeedEffVsSimHits",0,1);

  //  TCanvas *c1_06 = newCanvas("FigSTA/efficiencies/standAlone/sta_eff_vs_simhits","Muon reconstruction efficiences");  
  //  TH1F *h1_06 = adamDrawLoop(dir,"StandAloneMuonsUpdatedAtVtxEffVsSimHits",0,1);

  /////
  // tracker _2_
  /////
  TCanvas *c2_01 = newCanvas(dir+"FigGLB/efficiencies/tracker/tk_sim_eff_vs_eta","Muon reconstruction efficiences");
  TH1F *h2_01 = adamDrawLoop(dir,"GeneralTracksSimEffVsEta",0,1,1,1);

  TCanvas *c2_02 = newCanvas(dir+"FigGLB/efficiencies/tracker/tk_sim_eff_vs_phi","Muon reconstruction efficiences");
  TH1F *h2_02 = adamDrawLoop(dir,"GeneralTracksSimEffVsPhi",0,1,1,2);


  /////
  // tracker _22_
  /////
  
  TCanvas *c22_01 = newCanvas(dirTkmu+"FigGLB/efficiencies/trackerMuons/tkmu_sim_eff_vs_eta","Muon reconstruction efficiences");
  TH1F *h22_01 = adamDrawLoop(dirTkmu,"TrackerMuonsTrackerOnlySimEffVsEta",0,1,1,1);
  
  TCanvas *c22_02 = newCanvas(dirTkmu+"FigGLB/efficiencies/trackerMuons/tkmu_sim_eff_vs_phi","Muon reconstruction efficiences");
  TH1F *h22_02 = adamDrawLoop(dirTkmu,"TrackerMuonsTrackerOnlySimEffVsPhi",0,1,1,2);
  

  /////
  // GlobalMuons _3_
  /////

  TCanvas *c3_01 = newCanvas(dir+"FigGLB/efficiencies/global/glb_sim_eff_vs_eta","Muon reconstruction efficiences");
  TH1F *h3_01 = adamDrawLoop(dir,"GlobalMuonsSimEffVsEta",0,1,1,1);

  TCanvas *c3_02 = newCanvas(dir+"FigGLB/efficiencies/global/glb_sim_eff_vs_phi","Muon reconstruction efficiences");
  TH1F *h3_02 = adamDrawLoop(dir,"GlobalMuonsSimEffVsPhi",0,1,1,2);

  TCanvas *c3_03 = newCanvas(dir+"FigGLB/efficiencies/global/glb_sta_eff_vs_eta","Muon reconstruction efficiences");
  TH1F *h3_03 = adamDrawLoop(dir,"GlobalMuonsStandAloneMuonsUpdatedAtVtxEffVsEta",0,1,1,1);

  TCanvas *c3_04 = newCanvas(dir+"FigGLB/efficiencies/global/glb_sta_eff_vs_phi","Muon reconstruction efficiences");
  TH1F *h3_04 = adamDrawLoop(dir,"GlobalMuonsStandAloneMuonsUpdatedAtVtxEffVsPhi",0,1,1,2);

  TCanvas *c3_05 = newCanvas(dir+"FigGLB/efficiencies/global/glb_tk_eff_vs_eta","Muon reconstruction efficiences");
  TH1F *h3_05 = adamDrawLoop(dir,"GlobalMuonsGeneralTracksEffVsEta",0,1,1,1);

  TCanvas *c3_06 = newCanvas(dir+"FigGLB/efficiencies/global/glb_tk_eff_vs_phi","Muon reconstruction efficiences");
  TH1F *h3_06 = adamDrawLoop(dir,"GlobalMuonsGeneralTracksEffVsPhi",0,1,1,2);

  //////
  // Unbinned
  /////
  TCanvas *c4_01 = newCanvas(dir+"FigGLB/efficiencies/eff_vs_eta_0_500","Muon reconstruction efficiences");
  TH1F *h4_01 = adamDrawLoop(dir,"SimEffVsEta",2,1,1,1);

  TCanvas *c4_02 = newCanvas(dir+"FigGLB/efficiencies/eff_vs_phi_0_500","Muon reconstruction efficiences");
  TH1F *h4_02 = adamDrawLoop(dir,"SimEffVsPhi",2,1,1,2);

  TCanvas *c4_03 = newCanvas(dir+"FigGLB/efficiencies/eff_vs_pt_0_500","Muon reconstruction efficiences");
  TH1F *h4_03 = adamDrawLoop(dir,"SimEffVsPt",2,1,1,3);

  TCanvas *c4_04 = newCanvas(dir+"FigGLB/efficiencies/eff_vs_pt_Zoom_0_500","Muon reconstruction efficiences");
  TH1F *h4_04 = adamDrawLoop(dir,"SimEffVsPt",2,1,1,3);
  h4_04->GetXaxis()->SetRangeUser(0.,20.);

  //printCanvasesType(".eps");
  //printCanvasesType(".pdf");
  //printCanvasesPS("allPostMuEff.ps");
  
}
