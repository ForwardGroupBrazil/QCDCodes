void TriggerEfficiency()
{
  //TFile *inf = TFile::Open("data/flatTree_Jet.root");
  //TTree *tr  = (TTree*)inf->Get("Hbb/events"); 
  TChain *ch = new TChain("chain","chain");
  ch->Add("root://eoscms//eos/cms/store/cmst3/user/kkousour/VBF-CMG/flatTree_Jet-Run2012A.root/Hbb/events");
  ch->Add("root://eoscms//eos/cms/store/cmst3/user/kkousour/VBF-CMG/flatTree_JetMon-Run2012B.root/Hbb/events");
  ch->Add("root://eoscms//eos/cms/store/cmst3/user/kkousour/VBF-CMG/flatTree_JetMon-Run2012C.root/Hbb/events");
  ch->Add("root://eoscms//eos/cms/store/cmst3/user/kkousour/VBF-CMG/flatTree_JetMon-Run2012D.root/Hbb/events");
  const int N = 9;
  TH1F *h[N],*hREF[N];
  TString OBS[N]  = {"dEtaqq","jetPt[0]","jetPt[1]","jetPt[2]","jetPt[3]","mqq","jetBtag[btagIdx[0]]","jetBtag[btagIdx[0]]","mbb"};
  TString XTITLE[N] = {"|#Delta#eta_{qq}|", "Leading jet p_{T} (GeV)","Second jet p_{T} (GeV)",
                       "Third jet p_{T} (GeV)", "Fourth jet p_{T} (GeV)","mqq (GeV)","CSV[0]","CSV[1]","mbb (GeV)"};
  int REF_TRIG[N] = {9,9,9,9,9,9,9,9,9};
  int NBINS[N]    = {10,12,12,9,16,5,10,10,5};
  double XMIN[N]  = {0,20,20,40,20,250,0,0,50};
  double XMAX[N]  = {6,140,140,130,100,500,1.001,1.001,200};
  double CUT[N]   = {3.0,90,70,60,40,300,-1,-1,0};
  char name[1000];
  TCanvas *can[N];
  TEfficiency *hEff[N];
  TString SEL,ss;
  for(int i=0;i<9;i++) {
    h[i] = new TH1F(OBS[i],OBS[i],NBINS[i],XMIN[i],XMAX[i]);
    hREF[i] = new TH1F("REF_"+OBS[i],"REF_"+OBS[i],NBINS[i],XMIN[i],XMAX[i]);
    sprintf(name,"triggerResult[%d]==1",REF_TRIG[i]);
    TCut cutREF(name);
    sprintf(name,"triggerResult[%d]==1 && (triggerResult[0]==1 || triggerResult[1]==1)",REF_TRIG[i]);
    TCut cutSIG(name);
    if (i==0) {
      SEL = "jetPt[0]>90 && jetPt[1]>70 && jetPt[2]>60 && jetPt[3]>40 && mqq>300";
    }
    else if (i==1) {
      SEL = "dEtaqq>3 && jetPt[1]>70 && jetPt[2]>60 && jetPt[3]>40 && mqq>300";
    }
    else if (i==2) {
      SEL = "dEtaqq>3 && jetPt[0]>90 && jetPt[2]>60 && jetPt[3]>40 && mqq>300";
    }  
    else if (i==3) {
      SEL = "dEtaqq>3 && jetPt[0]>90 && jetPt[1]>70 && jetPt[3]>40 && mqq>300";
    }
    else if (i==4) {
      SEL = "dEtaqq>3 && jetPt[0]>90 && jetPt[1]>70 && jetPt[2]>60 && mqq>300";
    }
    else if (i==5) {
      SEL = "dEtaqq>3 && jetPt[0]>90 && jetPt[1]>70 && jetPt[2]>60 && jetPt[3]>40";
    }
    else if (i==6 || i==7) {
      SEL = "dEtaqq>3 && jetPt[0]>90 && jetPt[1]>70 && jetPt[2]>60 && jetPt[3]>40 && mqq>300";
    }
    else {
      SEL = "dEtaqq>3 && jetPt[0]>90 && jetPt[1]>70 && jetPt[2]>60 && jetPt[3]>40 && jetBtag[btagIdx[0]]>0.898 && jetBtag[btagIdx[1]]>0.898 && mqq>300";
    }
    TCut cutSEL(SEL);
    can[i] = new TCanvas("Eff_"+OBS[i],"Eff_"+OBS[i],900,600);
    gPad->SetGridx();
    gPad->SetGridy();
    ss = OBS[i]+">>"+"REF_"+OBS[i];
    ch->Draw(ss,cutREF+cutSEL);
    ss = OBS[i]+">>"+OBS[i];
    ch->Draw(OBS[i]+">>"+OBS[i],cutSIG+cutSEL); 
    hEff[i] = new TEfficiency(("Efficiency_"+OBS[i]).Data(),("Efficiency_"+OBS[i]+";"+XTITLE[i]+";Trigger Efficiency").Data(),NBINS[i],XMIN[i],XMAX[i]);
    hEff[i]->SetPassedHistogram(*h[i],"f");
    hEff[i]->SetTotalHistogram(*hREF[i],"f");
    hEff[i]->Draw();
    gPad->Update();
    if (OBS[i] != "mbb") {
      TLine *ln = new TLine(CUT[i],gPad->GetFrame()->GetY1(),CUT[i],gPad->GetFrame()->GetY2());
      ln->SetLineColor(kRed);
      ln->SetLineWidth(2);
      ln->SetLineStyle(2);
      ln->Draw("same");
    }
    else {
      TF1 *fit = new TF1("fit","[0]",gPad->GetFrame()->GetX1(),gPad->GetFrame()->GetX2());
      fit->SetParameter(0,0.5);
      fit->SetLineColor(kRed);
      fit->SetLineWidth(2);
      fit->SetLineStyle(2);
      hEff[i]->Fit(fit,"RQ+"); 
      fit->Draw("same");
    }
  }
}
