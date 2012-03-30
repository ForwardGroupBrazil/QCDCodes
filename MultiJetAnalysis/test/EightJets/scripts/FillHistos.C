void FillHistos(TString FileName, bool UseWeights)
{
  TFile *inf,*outf;
  TTree *tr;
  cout<<"Reading File: "<<FileName<<endl;  
  inf = TFile::Open(FileName);    
  outf = new TFile("Histo_"+FileName,"RECREATE");
  //----- booking histograms ---------------------- 
  TH1F *hBDT   = new TH1F("BDT","BDT",100,-1,1); 
  TH1F *hMLP   = new TH1F("MLP","MLP",30,-0.2,1.1);
  TH1F *hCos   = new TH1F("Cos","Cos",100,-1,1);  
  TH1F *hHT    = new TH1F("HT","HT",70,0,3.5);
  TH1F *hM8J   = new TH1F("M8J","M8J",70,0,7);
  TH1F *hM4J   = new TH1F("M4J","M4J",50,0,2);
  TH1F *hHT4J  = new TH1F("HT4J","HT4J",100,0,2.0);  
  TH1F *hPT4J  = new TH1F("PT4J","PT4J",60,0,1.5);
  TH1F *hEta4J = new TH1F("Eta4J","Eta4J",50,-5,5);  
  TH1F *hM2J   = new TH1F("M2J","M2J",200,0,1.0);
  TH1F *hB4J   = new TH1F("B4J","B4J",100,0,2);
  TH1F *hB2J   = new TH1F("B2J","B2J",100,0,2); 
  TH1F *hX[8],*hPt[8],*hEta[8],*hPhi[8],*hBeta[8];
  TH1F *hBDTCut = new TH1F("BDTCut","BDTCut",100,-1,1); 
  TH1F *hMLPCut = new TH1F("MLPCut","MLPCut",40,-0.2,1.2);
  TH1F *hCosCut   = new TH1F("CosCut","CosCut",100,-1,1);  
  TH1F *hHTCut    = new TH1F("HTCut","HTCut",70,0,3.5);
  TH1F *hM8JCut   = new TH1F("M8JCut","M8JCut",70,0,7);
  TH1F *hM4JCut   = new TH1F("M4JCut","M4JCut",50,0,2);
  TH1F *hHT4JCut  = new TH1F("HT4JCut","HT4JCut",100,0,2.0);  
  TH1F *hPT4JCut  = new TH1F("PT4JCut","PT4JCut",60,0,1.5);
  TH1F *hEta4JCut = new TH1F("Eta4JCut","Eta4JCut",50,-5,5);  
  TH1F *hM2JCut   = new TH1F("M2JCut","M2JCut",200,0,1.0);
  TH1F *hB4JCut   = new TH1F("B4JCut","B4JCut",100,0,2);
  TH1F *hB2JCut   = new TH1F("B2JCut","B2JCut",100,0,2); 
  TH1F *hXCut[8],*hPtCut[8],*hEtaCut[8],*hPhiCut[8],*hBetaCut[8];
  char name[1000];  
  for(int i=0;i<8;i++) {
    sprintf(name,"JetPtRatio%d",i);
    hX[i] = new TH1F(name,name,100,0,1.00001);  
    hX[i]->Sumw2();  
    sprintf(name,"JetPt%d",i);
    hPt[i] = new TH1F(name,name,100,0,1000);  
    hPt[i]->Sumw2();
    sprintf(name,"JetEta%d",i);
    hEta[i] = new TH1F(name,name,50,-5,5);  
    hEta[i]->Sumw2();
    sprintf(name,"JetPhi%d",i);
    hPhi[i] = new TH1F(name,name,100,-3.15,3.15);  
    hPhi[i]->Sumw2();
    sprintf(name,"JetBeta%d",i);
    hBeta[i] = new TH1F(name,name,100,0,1.0001);  
    hBeta[i]->Sumw2();
    sprintf(name,"JetPtRatioCut%d",i);
    hXCut[i] = new TH1F(name,name,100,0,1.00001);  
    hXCut[i]->Sumw2();  
    sprintf(name,"JetPtCut%d",i);
    hPtCut[i] = new TH1F(name,name,100,0,1000);  
    hPtCut[i]->Sumw2();
    sprintf(name,"JetEtaCut%d",i);
    hEtaCut[i] = new TH1F(name,name,50,-5,5);  
    hEtaCut[i]->Sumw2();
    sprintf(name,"JetPhiCut%d",i);
    hPhiCut[i] = new TH1F(name,name,100,-3.15,3.15);  
    hPhiCut[i]->Sumw2();
    sprintf(name,"JetBetaCut%d",i);
    hBetaCut[i] = new TH1F(name,name,100,0,1.0001);  
    hBetaCut[i]->Sumw2();
  }
  hBDT->Sumw2();
  hMLP->Sumw2();
  hCos->Sumw2();  
  hHT->Sumw2();
  hM8J->Sumw2();
  hM4J->Sumw2();  
  hHT4J->Sumw2();
  hPT4J->Sumw2();
  hEta4J->Sumw2();  
  hM2J->Sumw2();
  hB4J->Sumw2();
  hB2J->Sumw2();
  hBDTCut->Sumw2();
  hMLPCut->Sumw2();
  hCosCut->Sumw2();  
  hHTCut->Sumw2();
  hM8JCut->Sumw2();
  hM4JCut->Sumw2();  
  hHT4JCut->Sumw2();
  hPT4JCut->Sumw2();
  hEta4JCut->Sumw2();  
  hM2JCut->Sumw2();
  hB4JCut->Sumw2();
  hB2JCut->Sumw2();
  //----- tree variables ---------------------- 
  float BDT,MLP,m2jAve,m2jSig,m4jAve,m4jBalance,m8j,ht,cosThetaStar,wt(1.0);
  float ht4j[2],pt4j[2],eta4j[2];
  float pt[8],eta[8],phi[8],beta[8];
  //----- tree branches -----------------------
  tr = (TTree*)inf->Get("multijets/events");  
  tr->SetBranchAddress("BDT",&BDT);
  tr->SetBranchAddress("MLP",&MLP);
  tr->SetBranchAddress("ht",&ht);
  tr->SetBranchAddress("m8j",&m8j); 
  tr->SetBranchAddress("m4jAve",&m4jAve);
  tr->SetBranchAddress("m4jBalance",&m4jBalance);  
  tr->SetBranchAddress("m2jAve",&m2jAve);
  tr->SetBranchAddress("m2jSig",&m2jSig);
  tr->SetBranchAddress("cosThetaStar",&cosThetaStar);
  tr->SetBranchAddress("ht4j",&ht4j);
  tr->SetBranchAddress("pt4j",&pt4j);
  tr->SetBranchAddress("eta4j",&eta4j);
  tr->SetBranchAddress("pt",&pt);
  tr->SetBranchAddress("eta",&eta);
  tr->SetBranchAddress("phi",&phi);
  tr->SetBranchAddress("beta",&beta);
  if (UseWeights) {  
    tr->SetBranchAddress("wt",&wt); 
  }    
  //----- loop over tree entries ----------------
  for(int iev=0;iev<tr->GetEntries();iev++) {
    tr->GetEntry(iev);
    bool preselection = (ht > 850 && pt[7] > 30);
    if (!preselection) continue;
    hBDT->Fill(BDT,wt);
    hMLP->Fill(MLP,wt);
    hHT->Fill(ht/1000,wt);
    hM8J->Fill(m8j/1000,wt);
    hM4J->Fill(m4jAve/1000,wt);
    hM2J->Fill(m2jAve/1000,wt);
    hHT4J->Fill(ht4j[0]/1000,wt);
    hHT4J->Fill(ht4j[1]/1000,wt);  
    hPT4J->Fill(pt4j[0]/1000,wt);
    hPT4J->Fill(pt4j[1]/1000,wt);
    hEta4J->Fill(eta4j[0],wt);
    hEta4J->Fill(eta4j[1],wt);  
    hCos->Fill(cosThetaStar,wt);
    hB4J->Fill(m4jBalance,wt);
    hB2J->Fill(m2jSig/m2jAve,wt);
    for(int j=0;j<8;j++) {
      hPt[j]->Fill(pt[j],wt);
      hEta[j]->Fill(eta[j],wt);
      hPhi[j]->Fill(phi[j],wt);
      hBeta[j]->Fill(beta[j],wt);
      hX[j]->Fill(pt[j]/pt[0],wt);  
    }// jet loop
    if (MLP > 0.6) {
      hBDTCut->Fill(BDT,wt);
      hMLPCut->Fill(MLP,wt);
      hHTCut->Fill(ht/1000,wt);
      hM8JCut->Fill(m8j/1000,wt);
      hM4JCut->Fill(m4jAve/1000,wt);
      hM2JCut->Fill(m2jAve/1000,wt);
      hHT4JCut->Fill(ht4j[0]/1000,wt);
      hHT4JCut->Fill(ht4j[1]/1000,wt);  
      hPT4JCut->Fill(pt4j[0]/1000,wt);
      hPT4JCut->Fill(pt4j[1]/1000,wt);
      hEta4JCut->Fill(eta4j[0],wt);
      hEta4JCut->Fill(eta4j[1],wt);  
      hCosCut->Fill(cosThetaStar,wt);
      hB4JCut->Fill(m4jBalance,wt);
      hB2JCut->Fill(m2jSig/m2jAve,wt);
      for(int j=0;j<8;j++) {
        hPtCut[j]->Fill(pt[j],wt);
        hEtaCut[j]->Fill(eta[j],wt);
        hPhiCut[j]->Fill(phi[j],wt);
        hBetaCut[j]->Fill(beta[j],wt);
        hXCut[j]->Fill(pt[j]/pt[0],wt);  
      }// jet loop
    }
  }// tree loop
  outf->Write();  
}  