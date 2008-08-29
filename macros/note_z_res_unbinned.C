#include <vector>
#include <TString>

void note_z_res_unbinned(TString &dir="./"){
  gROOT->LoadMacro("adamStyles.C");
  gROOT->LoadMacro("adamCanvasMacros.C");
  gROOT->LoadMacro("adamDrawMacros.C");
  gROOT->LoadMacro("adamDrawLoop.C");

  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();

  gStyle->SetTitleFillColor(0);

  gStyle->SetOptLogy(1);

  // gLegend = new TLegend(.4, .85, 0.6, .95, "");

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
			"res_qpt",
			"pull_d0",
                        "pull_dsz"};
  TString yTitle[] = {"#sigma(#delta d_{xy})[cm]",
		      "#sigma(#delta d_{z}) [cm]",
		      "#sigma(#delta #eta)",
		      "#sigma(#delta #phi)[rad]",
		      "#sigma((q/p_{t})/(q/p_{t}))",
		      "#sigma(#delta d_{0}) [cm]",
		      "#sigma(#delta d_{sz}) [cm]"};
  TString histo[] = {"ResoDxy",
		     "ResoDz",
		     "ResoEta",
		     "ResoPhi",
		     "ResoQOverPt",
		     "ResoD0",
		     "ResoDsz"};
  TString xAxis[]  = {"eta","phi","pt"}; 
  TString XAxis[]  = {"Eta","Phi","Pt"}; 
  
  int canvasCounter = 01;
  TCanvas * canvas;
  TH1 * h;

  for (int iLevel=0;iLevel < 4; ++iLevel) {
    for (int iQuantity=0;iQuantity < 7; ++iQuantity) {
      for (int iXaxis=0;iXaxis < 3; ++iXaxis) {
	canvas = newCanvas(dir+"FigParam/resolutions/"+level[iLevel]+"_"+quantity[iQuantity]+"_vs_"+xAxis[iXaxis]+"_unbinned","Resolutions "+levelName[iLevel]+" "+histo[iQuantity]);
	int myOpt = (iXaxis==2) ? 0 : 1;
	h = adamDrawLoop(dir,levelName[iLevel]+histo[iQuantity]+"Vs"+XAxis[iXaxis],999,myOpt,2,iXaxis+1);
	h->GetYaxis()->SetTitle(yTitle[iQuantity]);
	if(iXaxis==2) h->Rebin(10);
      }
    }
  }
  
  printCanvasesType(".eps");
  printCanvasesType(".pdf");
  printCanvasesPS(dir+"allPostMuResUnbinned.ps");

}
