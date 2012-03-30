using namespace RooFit;
void Fit()
{
  gROOT->ForceStyle();
  TFile *fDat = TFile::Open("Histo_flatTree_data_tmva.root");
  TFile *fBkg = TFile::Open("Histo_flatTree_qcd_weights_tmva.root");
  TFile *fSgn = TFile::Open("Histo_flatTree_Col800Sig200_weights_tmva.root");
  
  TH1 *hDat = (TH1*)fDat->Get("MLP");
  TH1 *hBkg = (TH1*)fBkg->Get("MLP");
  TH1 *hSgn = (TH1*)fSgn->Get("MLP");
  double k_factor = 1.3;
  double lumi = 5000;
  hBkg->Scale(k_factor*lumi);
  hSgn->Scale(lumi);
  
  RooRealVar x("x","x",-0.2,1.2);
  
  RooDataHist data("data","dataset with x",x,hDat);
  RooDataHist bkg("qcd","bkg with x",x,hBkg);
  RooDataHist sgn("signal","sgn with x",x,hSgn);
  
  RooHistPdf bkgPDF("bkgPDF","bkgPDF",x,bkg,0);
  RooHistPdf sgnPDF("sgnPDF","sgnPDF",x,sgn,0);
  
  RooRealVar f("f","f",0,0.,1.);
  
  RooAddPdf model("model","model",RooArgList(sgnPDF,bkgPDF),RooArgList(f));
  
  RooFitResult* r = model.fitTo(data,Save());
  r->Print("v");
  double N = hDat->GetEntries();
  double B = hBkg->Integral();
  double S = hSgn->Integral();
  cout<<N<<" "<<B<<" "<<S<<endl;
  
  RooPlot* frame1 = x.frame();  
  data.plotOn(frame1);
  model.plotOn(frame1,Components("sgnPDF*"),LineStyle(kDashed));
  model.plotOn(frame1);
  
  //cout<<frame1->chiSquare()<<endl;
  RooHist* hresid = frame1->residHist();
  RooHist* hpull = frame1->pullHist();
  RooPlot* frame2 = x.frame(Title("Residual Distribution"));
  frame2->addPlotable(hresid,"P") ;
  
  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* frame3 = x.frame(Title("Pull Distribution"));
  frame3->addPlotable(hpull,"P");
  
  TCanvas* c = new TCanvas("rf109_chi2residpull","rf109_chi2residpull",800,600);
  c->Divide(1,2);
  c->cd(1); 
  gPad->SetLogy();
  frame1->SetMaximum(1e+5);
  frame1->SetMinimum(1e-2);
  frame1->Draw();

  c->cd(2); 
  frame3->Draw();
  
  TFile *out = new TFile("eightJet_input.root","RECREATE");
  out->cd();
  hDat->Write("data_obs");
  hBkg->Write("background");
  hSgn->Write("signal");
}

