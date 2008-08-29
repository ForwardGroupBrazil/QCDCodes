#include <vector>
#include <TString>

void note_z_pull_unbinned(TString &dir="./"){
  gROOT->LoadMacro("adamStyles.C");
  gROOT->LoadMacro("adamCanvasMacros.C");
  gROOT->LoadMacro("adamDrawMacros.C");
  gROOT->LoadMacro("adamDrawLoop.C");

  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();

  gStyle->SetTitleFillColor(0);

  gStyle->SetOptLogy(0);

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
  TString quantity2[] = {"Pull_dxy",
			 "Pull_dz",
			 "Pull_eta",
			 "Pull_phi",
			 "Pull_qpt",
			 "Pull_d0",
			 "Pull_dsz"};
  TString yTitle[] = {"#sigma_{pull}(d_{xy})",
		      "#sigma_{pull}(d_{z})",
		      "#sigma_{pull}(#eta)",
		      "#sigma_{pull}(#phi)",
		      "#sigma_{pull}((q/p_{t})/(q/p_{t}))",
		      "#sigma_{pull}(d_{0})",
		      "#sigma_{pull}(d_{sz})"};
  TString yTitle2[] = {"Pull(d_{xy})",
		       "Pull(d_{z})",
		       "Pull(#eta)",
		       "Pull(#phi)",
		       "Pull((q/p_{t})/(q/p_{t}))",
		       "Pull(d_{0})",
		       "Pull(d_{sz})"};
  TString histo[] = {"PullDxy",
		     "PullDz",
		     "PullEta",
		     "PullPhi",
		     "PullQOverPt",
		     "PullD0",
		     "PullDsz"};

  TString xAxis[]  = {"eta","phi","pt"}; 
  TString XAxis[]  = {"Eta","Phi","Pt"}; 
  
  int canvasCounter = 01;
  TCanvas * canvas;
  TH1 * h;

  /////
  // For the note
  /////

  for (int iLevel=0;iLevel < 4; ++iLevel) {
    for (int iQuantity=0;iQuantity < 7; ++iQuantity) {
	canvas = newCanvas(dir+"FigParam/pulls/"+level[iLevel]+"_"+quantity2[iQuantity]+"_unbinned","Pulls "+levelName[iLevel]+" "+quantity2[iQuantity]);
	h = adamDrawLoop(dir,levelName[iLevel]+histo[iQuantity],999,0,1,0);
	h->GetYaxis()->SetTitle(yTitle2[iQuantity]);
	//if(histo[iQuantity] == "PullEta") h->Rebin(10);
	h->Fit("gaus");
	TF1 *myfunc = h->GetFunction("gaus");
	h->SetLineWidth(2);
	myfunc->SetLineWidth(2);
    }
  }

  /////
  // For completeness
  /////
  /*
  for (int iLevel=0;iLevel < 4; ++iLevel) {
    for (int iQuantity=0;iQuantity < 7; ++iQuantity) {
      for (int iXaxis=0;iXaxis < 3; ++iXaxis) {	
	canvas = newCanvas(dir+"FigParam/pulls/"+level[iLevel]+"_"+quantity[iQuantity]+"_vs_"+xAxis[iXaxis],"Pulls "+levelName[iLevel]+" "+quantity[iQuantity]);
	h = adamDrawLoop(dir,levelName[iLevel]+histo[iQuantity]+"Vs"+XAxis[iXaxis],999,1,2,iXaxis+1);
	h->GetYaxis()->SetTitle(yTitle[iQuantity]);
      }
    }
  }
  */


  //printCanvasesType(".eps");
  //printCanvasesType(".pdf");
  //printCanvasesPS(dir+"allPostMuPullUnbinned.ps");
  
}
