#include <vector>
#include <TString>

void note_z_res_inclusive(TString &dir="./"){
  gROOT->LoadMacro("adamStyles.C");
  gROOT->LoadMacro("adamCanvasMacros.C");
  gROOT->LoadMacro("adamDrawMacros.C");
  gROOT->LoadMacro("adamDrawLoop.C");

  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();

  gStyle->SetTitleFillColor(0);

  gStyle->SetOptLogy(1);

  dir = "/scratch/scratch96/a/aeverett/reReco_2011a/postMu/";
  TString dirTkmu = "/scratch/scratch96/a/aeverett/reReco_2011a/postTkmu/";

  TString level[] = {"muonSeed/seed",
		     "standAlone/sta",
		     "global/glb",
		     "tracker/tk"};
  TString levelName[] = {"MuonSeed",
			 "StandAloneMuonsUpdatedAtVtx",
			 "GlobalMuons", 
			 "GeneralTracks"};
  TString quantity[] = {"res_dxy",
			"res_dz",
			"res_eta",
			"res_phi",
			"res_qpt"};
  TString yTitle[] = {"#sigma(#delta d_{xy})[cm]",
		      "#sigma(#delta d_{z}) [cm]",
		      "#sigma(#delta #eta)",
		      "#sigma(#delta #phi)[rad]",
		      "#sigma((q/p_{t})/(q/p_{t}))"};
  TString histo[] = {"ResoDxy",
		     "ResoDz",
		     "ResoEta",
		     "ResoPhi",
		     "ResoQOverPt"};
  TString xAxis[]  = {"eta","phi"}; 
  TString XAxis[]  = {"Eta","Phi"}; 
  
  int canvasCounter = 01;
  TCanvas * canvas;
  TH1 * h;
  for (int iLevel=0;iLevel < 4; ++iLevel) {
    for (int iQuantity=0;iQuantity < 5; ++iQuantity) {
      for (int iXaxis=0;iXaxis < 2; ++iXaxis) {
	canvas = newCanvas(dir+"FigParam/resolutions/"+level[iLevel]+"_"+quantity[iQuantity]+"_vs_"+xAxis[iXaxis],"Resolutions "+levelName[iLevel]);
	h = adamDrawLoop(dir,levelName[iLevel]+histo[iQuantity]+"Vs"+XAxis[iXaxis],0,1,2,iXaxis+1);
	h->GetYaxis()->SetTitle(yTitle[iQuantity]);
	h->SetMinimum(0.01);
	h->SetMaximum(100.0);
      }
    }
  }


  printCanvasesType(".eps");
  //printCanvasesType(".pdf");
  printCanvasesPS("allPostMuRes.ps");

}
