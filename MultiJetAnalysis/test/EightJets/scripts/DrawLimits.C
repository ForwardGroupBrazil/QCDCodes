void DrawLimits()
{
  gROOT->ForceStyle();
  const int N = 11; 
  double signal[N][6];
  double xsecFull[N] = {0.13660E+02,0.31137E+01,0.81790E+00,0.23289E+00,0.71177E-01,0.22823E-01,0.75514E-02,0.25646E-02,0.88668E-03,0.30775E-03,0.10714E-03};
  double xsec[N] = {4.755,1.774,0.5974,0.1903,0.0618,0.02062,0.00695605,0.002393,0.0008348,0.0002932,0.0001026};
  double bayes[N] = {0.0193,0.0769,0.0811,0.0848,0.1316,0.0771327,0.245052,0.4345,1.15834,2.00986,6.989};
  double mass[N] = {500,600,700,800,900,1000,1100,1200,1300,1400,1500};
  double Nevt[N] = {25000,50000,25000,25000,50000,25000,25000,50000,25000,25000,25000};
  TGraphAsymmErrors *g1,*g2,*gObs,*gExp;
  double x[N],y1[N],y2[N],yExp[N],yObs[N],yBayes[N],yModel[N],yAcc[N];
  double exl[N],ey1l[N],ey2l[N],eyExpl[N],eyObsl[N],eyAcc[N];
  double exh[N],ey1h[N],ey2h[N],eyExph[N],eyObsh[N];
  double limit;
  char name[1000];
  for(int i=0;i<N;i++) {
    sprintf(name,"limits/higgsCombineColoronMLP.Asymptotic.mH%1.0f.root",mass[i]);
    TFile *inf = TFile::Open(name);
    TTree *tr = (TTree*)inf->Get("limit");
    tr->SetBranchAddress("limit",&limit);
    for(int iev=0;iev<tr->GetEntries();iev++) {
      tr->GetEntry(iev);
      signal[i][iev] = limit;
    }
    x[i] = mass[i];
    exl[i] = 0;
    exh[i] = 0;
    yAcc[i] = xsec[i]/xsecFull[i];
    eyAcc[i] = sqrt(yAcc[i]*(1-yAcc[i])/Nevt[i]);
    y1[i] = signal[i][2]*xsec[i];
    y2[i] = signal[i][2]*xsec[i];
    yExp[i] = signal[i][2]*xsec[i];
    yObs[i] = signal[i][5]*xsec[i];
    ey1l[i] = (signal[i][2]-signal[i][1])*xsec[i];
    ey2l[i] = (signal[i][2]-signal[i][0])*xsec[i];
    eyExpl[i] = 0;
    eyObsl[i] = 0;
    ey1h[i] = (signal[i][3]-signal[i][2])*xsec[i];
    ey2h[i] = (signal[i][4]-signal[i][2])*xsec[i];
    eyExph[i] = 0;
    eyObsh[i] = 0;
    yBayes[i] = bayes[i]*xsec[i];
    yModel[i] = xsec[i];
  }
  TGraph *gBayes = new TGraph(N,mass,yBayes);
  TGraph *gModel = new TGraph(N,mass,yModel);
  TGraphErrors *gAcc   = new TGraphErrors(N,mass,yAcc,exl,eyAcc);
  g1 = new TGraphAsymmErrors(N,x,y1,exl,exh,ey1l,ey1h);
  g2 = new TGraphAsymmErrors(N,x,y2,exl,exh,ey2l,ey2h);
  gObs = new TGraphAsymmErrors(N,x,yObs,exl,exh,eyObsl,eyObsh);
  gExp = new TGraphAsymmErrors(N,x,yExp,exl,exh,eyExpl,eyExph);
  
  g2->SetFillColor(kYellow);
  g2->SetLineColor(kYellow);
  g1->SetFillColor(kGreen);
  g1->SetLineColor(kGreen);
  gExp->SetLineWidth(2);
  gExp->SetLineStyle(9);
  gObs->SetLineWidth(2);
  gObs->SetLineStyle(1);
  gObs->SetLineColor(kBlue+1);
  gObs->SetMarkerColor(kBlue+1);
  gObs->SetMarkerStyle(21);
  gObs->SetMarkerSize(1.5);
  gModel->SetLineStyle(2);
  gModel->SetLineWidth(2);
  gModel->SetLineColor(kRed);
  
  TCanvas *can = new TCanvas("Limits","Limits",900,600);
  gPad->SetLogy();
  gPad->SetGridy();
  gPad->SetGridx();
  g2->GetXaxis()->SetTitle("Coloron Mass (GeV)");
  g2->GetYaxis()->SetTitle("#sigma #times BR #times A (pb)");
  g2->SetMinimum(1e-4);
  g2->SetMaximum(5);
  g2->Draw("AE3");
  g1->Draw("sameE3");
  gExp->Draw("sameL");
  gObs->Draw("sameLP");
  gBayes->Draw("sameLP");
  gModel->Draw("sameL");
  TLegend *leg = new TLegend(0.5,0.65,0.9,0.9);
  leg->AddEntry(gObs,"Observed","LP");
  leg->AddEntry(gExp,"Expected Median","L");
  leg->AddEntry(g1,"Expected 1-#sigma","F");
  leg->AddEntry(g2,"Expected 2-#sigma","F");
  leg->AddEntry(gModel,"Coloron Model","L");
  leg->SetFillColor(0);
  //leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->Draw();
  
  TCanvas *canAcc = new TCanvas("Acceptance","Acceptance",900,600);
  gPad->SetGridx();
  gPad->SetGridy();
  gAcc->SetMinimum(0);
  gAcc->SetMaximum(1);
  gAcc->SetLineWidth(2);
  gAcc->Draw("ACP");
  gAcc->GetXaxis()->SetTitle("Coloron Mass (GeV)");
  gAcc->GetYaxis()->SetTitle("Acceptance");
}







