{
   TFile *_file0 = new TFile("DYM1020/allhistos16BinsStoyan1bbis/merge.root");
   TFile *_file1 = new TFile("DYM20/allhistos16BinsStoyan1bbis/merge.root");

   TString outputDir = "16BinsStoyan1a/";

   TFile *fileOut = new TFile("Stoyan1bbis.root","RECREATE");
   
   gROOT->SetStyle("Plain");
   gStyle->SetPalette(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   
   double w1 = 0.0243614;
   double w2 = 0.0121922;
   
   //TString dirs[] = {"genMuons","genMuonsBin1","genMuonsBin2","genMuonsBin3","genMuonsBin4","genMuonsBin5","genMuonsBin6","genMuonsBin7","genMuonsBin8","genMuonsBin9","genMuonsBin10","genMuonsBin11","genMuonsBin12","genMuonsBin13","genMuonsBin14","genMuonsBin15","genMuonsBin16"};//17
   TString dirs[] = {"genMuons"};//1
   TString histos[] = {"massGenZ_init","massGenDiMu_init","massGenZ_acc","massGenDiMu_acc","massGenZ_reco","massGenDiMu_reco"};//6
   
   for(int i=0; i<1; i++){
     cout << endl << "Looking in bin " + dirs[i] << " of " << outputDir <<endl;
     //temporary: in the analyzer, the REC and ACC histos were switched, so the indexes here are out of order

     TH1F* nInitGEN = (TH1F*)(_file0->Get(dirs[i]+"/"+histos[0])->Clone());
     TH1F* nInitFSR = (TH1F*)(_file0->Get(dirs[i]+"/"+histos[1])->Clone());
     TH1F* nAccGEN  = (TH1F*)(_file0->Get(dirs[i]+"/"+histos[2])->Clone());
     TH1F* nAccFSR  = (TH1F*)(_file0->Get(dirs[i]+"/"+histos[3])->Clone());
     TH1F* nRecoGEN = (TH1F*)(_file0->Get(dirs[i]+"/"+histos[4])->Clone());
     TH1F* nRecoFSR = (TH1F*)(_file0->Get(dirs[i]+"/"+histos[5])->Clone());
     
     TH1F* nInitGEN2 = (TH1F*)(_file1->Get(dirs[i]+"/"+histos[0])->Clone());
     TH1F* nInitFSR2 = (TH1F*)(_file1->Get(dirs[i]+"/"+histos[1])->Clone());
     TH1F* nAccGEN2  = (TH1F*)(_file1->Get(dirs[i]+"/"+histos[2])->Clone());
     TH1F* nAccFSR2  = (TH1F*)(_file1->Get(dirs[i]+"/"+histos[3])->Clone());
     TH1F* nRecoGEN2 = (TH1F*)(_file1->Get(dirs[i]+"/"+histos[4])->Clone());
     TH1F* nRecoFSR2 = (TH1F*)(_file1->Get(dirs[i]+"/"+histos[5])->Clone());
     
     nInitGEN->Add(nInitGEN, nInitGEN2, w1, w2);
     nInitFSR->Add(nInitFSR, nInitFSR2, w1, w2);
     nAccGEN ->Add(nAccGEN,  nAccGEN2,  w1, w2);
     nAccFSR ->Add(nAccFSR,  nAccFSR2,  w1, w2);
     nRecoGEN->Add(nRecoGEN, nRecoGEN2, w1, w2);
     nRecoFSR->Add(nRecoFSR, nRecoFSR2, w1, w2);
   
     TH1F* h_accGEN = (TH1F*) nInitGEN->Clone("h_accGEN");
     h_accGEN->SetTitle("Acceptance(GEN)");
     h_accGEN->Reset();
     h_accGEN->Sumw2();

     TH1F* h_accFSR = (TH1F*) nInitFSR->Clone("h_accFSR");
     h_accFSR->SetTitle("Acceptance(FSR)");
     h_accFSR->Reset();
     h_accFSR->Sumw2();

     TH1F* h_effGEN = (TH1F*) nInitGEN->Clone("h_effGEN");
     h_effGEN->SetTitle("Efficiency(GEN)");
     h_effGEN->Reset();
     h_effGEN->Sumw2();

     TH1F* h_effFSR = (TH1F*) nInitFSR->Clone("h_effFSR");
     h_effFSR->SetTitle("Efficiency(FSR)");
     h_effFSR->Reset();
     h_effFSR->Sumw2();

     TH1F* h_acceffGEN = (TH1F*) nInitGEN->Clone("h_acceffGEN");
     h_acceffGEN->SetTitle("Acceptance*Efficiency(GEN)");
     h_acceffGEN->Reset();
     h_acceffGEN->Sumw2();

     TH1F* h_acceffFSR = (TH1F*) nInitFSR->Clone("h_acceffFSR");
     h_acceffFSR->SetTitle("Acceptance*Efficiency(FSR)");
     h_acceffFSR->Reset();
     h_acceffFSR->Sumw2();
  
     int maxBins = nInitGEN->GetNbinsX();
     for(int bin=1; bin < maxBins+1; bin++) {
       cout << "Bin " << bin << " lower edge " << nInitGEN->GetBinLowEdge(bin)
	    << " upper edge " << nInitGEN->GetBinLowEdge(bin)+nInitGEN->GetBinWidth(bin) << endl;
       cout << "   nInitGEN " << nInitGEN->GetBinContent(bin) << endl;
       cout << "   nInitFSR " << nInitFSR->GetBinContent(bin) << endl;
       cout << "   nAccGEN  " << nAccGEN ->GetBinContent(bin) << endl;
       cout << "   nAccFSR  " << nAccFSR ->GetBinContent(bin) << endl;
       cout << "   nRecoGEN " << nRecoGEN->GetBinContent(bin) << endl;
       cout << "   nRecoFSR " << nRecoFSR->GetBinContent(bin) << endl;

       float accGEN = nInitGEN->GetBinContent(bin) ? ((float)(nAccGEN ->GetBinContent(bin)))/((float)(nInitGEN->GetBinContent(bin))) : 0.;
       float accGENerror = nInitGEN->GetBinContent(bin) && accGEN <= 1 ? sqrt(accGEN*(1-accGEN)/((float)(nInitGEN->GetBinContent(bin)))) : 0.;
       cout << "   Acc(GEN) " << accGEN << " +- " << accGENerror << endl;

       h_accGEN->SetBinContent(bin,accGEN);
       h_accGEN->SetBinError(bin,accGENerror);

       float accFSR = nInitFSR->GetBinContent(bin) ? ((float)(nAccFSR ->GetBinContent(bin)))/((float)(nInitFSR->GetBinContent(bin))) : 0.;
       float accFSRerror = nInitFSR->GetBinContent(bin) && accFSR <= 1 ? sqrt(accFSR*(1-accFSR)/((float)(nInitFSR->GetBinContent(bin)))) : 0.;
       cout << "   Acc(FSR) " << accFSR << " +- " << accFSRerror << endl;

       h_accFSR->SetBinContent(bin,accFSR);
       h_accFSR->SetBinError(bin,accFSRerror);

       float effGEN = nAccGEN->GetBinContent(bin) ? ((float)(nRecoGEN ->GetBinContent(bin)))/((float)(nAccGEN->GetBinContent(bin))) : 0.;
       float effGENerror = nAccGEN->GetBinContent(bin) && effGEN <= 1 ? sqrt(effGEN*(1-effGEN)/((float)(nAccGEN->GetBinContent(bin)))) : 0.;
       cout << "   Eff(GEN) " << effGEN << " +- " << effGENerror << endl;

       h_effGEN->SetBinContent(bin,effGEN);
       h_effGEN->SetBinError(bin,effGENerror);

       float effFSR = nAccFSR->GetBinContent(bin) ? ((float)(nRecoFSR ->GetBinContent(bin)))/((float)(nAccFSR->GetBinContent(bin))) : 0.;
       float effFSRerror = nAccFSR->GetBinContent(bin) && effFSR <= 1 ? sqrt(effFSR*(1-effFSR)/((float)(nAccFSR->GetBinContent(bin)))) : 0.;
       cout << "   Eff(FSR) " << effFSR << " +- " << effFSRerror << endl;

       h_effFSR->SetBinContent(bin,effFSR);
       h_effFSR->SetBinError(bin,effFSRerror);

       float acceffGEN = nInitGEN->GetBinContent(bin) ? ((float)(nRecoGEN ->GetBinContent(bin)))/((float)(nInitGEN->GetBinContent(bin))) : 0.;
       float acceffGENerror = nInitGEN->GetBinContent(bin) && acceffGEN <= 1 ? sqrt(acceffGEN*(1-acceffGEN)/((float)(nInitGEN->GetBinContent(bin)))) : 0.; 
       cout << "   Acc*Eff(GEN) " << acceffGEN << " +- " << acceffGENerror << endl;

       h_acceffGEN->SetBinContent(bin,acceffGEN);
       h_acceffGEN->SetBinError(bin,acceffGENerror);

       float acceffFSR = nInitFSR->GetBinContent(bin) ? ((float)(nRecoFSR ->GetBinContent(bin)))/((float)(nInitFSR->GetBinContent(bin))) : 0.;
       float acceffFSRerror = nInitFSR->GetBinContent(bin) && acceffFSR <= 1 ? sqrt(acceffFSR*(1-acceffFSR)/((float)(nInitFSR->GetBinContent(bin)))) : 0.; 
       cout << "   Acc*Eff(FSR) " << acceffFSR << " +- " << acceffFSRerror << endl;

       h_acceffFSR->SetBinContent(bin,acceffFSR);
       h_acceffFSR->SetBinError(bin,acceffFSRerror);

     }
     //h_accGEN->Draw();
     //h_accFSR->Draw();
   }
   fileOut->Write();
   fileOut->Close();
}
