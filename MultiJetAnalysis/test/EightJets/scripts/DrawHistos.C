void DrawHistos(TString NAME,int REBIN,bool LOG,TString XTITLE,double XMIN,double XMAX,bool PRINT,TString SIGNAL,TString TRAIN)
{
  gROOT->ForceStyle();
  float LUMI(4967);  
  TFile *fDAT = TFile::Open("Histo_flatTree_data_tmva"+SIGNAL+".root");
  TFile *fQCD = TFile::Open("Histo_flatTree_qcd_weights_tmva"+SIGNAL+".root");
  TFile *fSGN = TFile::Open("Histo_flatTree_"+SIGNAL+"_weights_tmva"+TRAIN+".root");
  
  TH1F *hDAT = (TH1F*)fDAT->Get(NAME);
  TH1F *hQCD = (TH1F*)fQCD->Get(NAME);
  TH1F *hSGN = (TH1F*)fSGN->Get(NAME);
  hDAT->Rebin(REBIN);
  hQCD->Rebin(REBIN);
  hSGN->Rebin(REBIN);
  TH1F *hDATRatio = (TH1F*)hDAT->Clone(NAME+"Ratio");
  
  hQCD->Scale(LUMI);
  hSGN->Scale(LUMI);
  TH1F *hSGNRatio = (TH1F*)hSGN->Clone(NAME+"Ratio");
  
  float norm = hDAT->Integral()/hQCD->Integral();
  hQCD->Scale(norm);
  cout<<"QCD k-factor: "<<norm<<endl;
  cout<<hQCD->Integral()<<endl;
  cout<<hSGN->Integral()<<endl;
  TH1F *hTOT = (TH1F*)hQCD->Clone(NAME+"TOT");
  hTOT->Add(hSGN);
  TH1F *hTOTRatio = (TH1F*)hTOT->Clone(NAME+"TOTRatio");
  hDATRatio->Divide(hQCD);
  hTOTRatio->Divide(hQCD);
  
  hQCD->SetFillColor(kYellow-10);
  hSGN->SetFillColor(kBlue-6);
  TH1F *hSGN2 = (TH1F*)hSGN->Clone("SGN2");
  hSGN2->SetLineWidth(2);
  hSGN2->SetLineColor(kBlue);
  hSGN2->SetFillColor(kWhite);
  
  hSGNRatio->SetLineColor(kBlue);
  hSGNRatio->SetLineWidth(2);
  
  THStack *hs = new THStack("hs","hs");
  hs->Add(hQCD);
  hs->Add(hSGN);

  
  
  TCanvas *can = new TCanvas("can_"+NAME+"_"+SIGNAL,"can_"+NAME+"_"+SIGNAL,900,600);
  if (LOG)
    gPad->SetLogy(); 
  hQCD->SetMinimum(0.5);
  hQCD->GetYaxis()->SetTitle("Events");
  hQCD->GetXaxis()->SetTitle(XTITLE);
  hQCD->GetXaxis()->SetRangeUser(XMIN,XMAX);
  hQCD->Draw("HIST");
  hDAT->Draw("same E");
  hSGN2->Draw("same HIST");
  TLegend *leg = new TLegend(0.65,0.6,0.9,0.9);
  leg->SetHeader("L = 5 fb^{-1}");
  leg->AddEntry(hDAT,"data","P");
  leg->AddEntry(hQCD,"QCD (#times 1.3)","F");
  leg->AddEntry(hSGN2,"signal","L");
  leg->SetBorderSize(0);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->Draw();
  
  TCanvas *canRatio = new TCanvas("can_"+NAME+"Ratio"+"_"+SIGNAL,"can_"+NAME+"Ratio"+"_"+SIGNAL,900,600);
  gPad->SetGridy();
  hDATRatio->GetYaxis()->SetTitle("Ratio to QCD");
  hDATRatio->GetXaxis()->SetTitle(XTITLE);
  hDATRatio->GetXaxis()->SetRangeUser(XMIN,XMAX);
  hDATRatio->SetMaximum(4.0);
  hDATRatio->SetMinimum(0.5);
  hTOTRatio->SetLineWidth(2);
  hTOTRatio->SetLineColor(kBlue);
  hDATRatio->Draw("E");
  hTOTRatio->Draw("same HIST ][");
  TLegend *leg1 = new TLegend(0.2,0.6,0.5,0.9);
  leg1->SetHeader("L = 5 fb^{-1}");
  leg1->AddEntry(hDATRatio,"data","P");
  leg1->AddEntry(hTOTRatio,"QCD+signal","L");
  leg1->SetBorderSize(0);
  leg1->SetLineColor(0);
  leg1->SetFillColor(0);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.05);
  leg1->Draw();
  
  TCanvas *canStack = new TCanvas("can_Stack"+NAME+"_"+SIGNAL,"can_Stack"+NAME+"_"+SIGNAL,900,600);
  if (LOG)
    gPad->SetLogy();
  hQCD->SetMinimum(0.5);
  hQCD->GetYaxis()->SetTitle("Events");
  hQCD->GetXaxis()->SetTitle(XTITLE);
  hQCD->GetXaxis()->SetRangeUser(XMIN,XMAX);
  hQCD->Draw("HIST");
  hs->Draw("same HIST");
  hDAT->Draw("same E");
  gPad->RedrawAxis();
  TLegend *leg = new TLegend(0.65,0.6,0.9,0.9);
  leg->SetHeader("L = 5 fb^{-1}");
  leg->AddEntry(hDAT,"data","P");
  leg->AddEntry(hQCD,"QCD (#times 1.3)","F");
  leg->AddEntry(hSGN,"signal","F");
  leg->SetBorderSize(0);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->Draw();
  
  if (PRINT) {
    can->Print("plots/"+TString(can->GetName())+".pdf");
    canRatio->Print("plots/"+TString(canRatio->GetName())+".pdf");
    canStack->Print("plots/"+TString(canStack->GetName())+".pdf");
  }
}