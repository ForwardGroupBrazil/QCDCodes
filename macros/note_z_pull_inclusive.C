#include <vector>
#include <TString>

void note_z_pull_inclusive(TString &dir="./"){
  gROOT->LoadMacro("adamStyles.C");
  gROOT->LoadMacro("adamCanvasMacros.C");
  gROOT->LoadMacro("adamDrawMacros.C");
  gROOT->LoadMacro("adamDrawLoop.C");

  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();

  gStyle->SetTitleFillColor(0);

  gStyle->SetOptLogy(0);

  //TLegend should be gLegend = new TLegend(.35, .13, 0.65, .33, "");

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
  
  TString quantity[] = {"pull_dxy",
			"pull_dz",
			"pull_eta",
			"pull_phi",
			"pull_qpt",
			"pull_d0",
			"pull_dsz"};
  TString yTitle[] = {"#sigma_{pull}(#delta d_{xy})",
		      "#sigma_{pull}(#delta d_{z})",
		      "#sigma_{pull}(#delta #eta)",
		      "#sigma_{pull}(#delta #phi)",
		      "#sigma_{pull}((q/p_{t})/(q/p_{t}))",
		      "#sigma_{pull}(d_{0})",
		      "#sigma_{pull}(d_{sz})"};
  TString histo[] = {"PullDxy",
		     "PullDz",
		     "PullEta",
		     "PullPhi",
		     "PullQOverPt",
		     "PullD0",
		     "PullDsz"};
  TString xAxis[]  = {"eta","phi"}; 
  TString XAxis[]  = {"Eta","Phi"}; 
  
  int canvasCounter = 01;
  TCanvas * canvas;
  TH1 * h;
  
  for (int iLevel=0;iLevel < 4; ++iLevel) {
    for (int iQuantity=0;iQuantity < 7; ++iQuantity) {
      for (int iXaxis=0;iXaxis < 2; ++iXaxis) {
	canvas = newCanvas(dir+"FigParam/pulls/"+level[iLevel]+"_"+quantity[iQuantity]+"_vs_"+xAxis[iXaxis],"Pulls");
	h = adamDrawLoop(dir,levelName[iLevel]+histo[iQuantity]+"Vs"+XAxis[iXaxis],0,1,2,iXaxis+1);
	h->GetYaxis()->SetTitle(yTitle[iQuantity]);
      }
    }
  }

  //printCanvasesType(".eps");
  //printCanvasesType(".pdf");
  printCanvasesPS(dir+"allPostMuPullInclusive.ps");
  
}
