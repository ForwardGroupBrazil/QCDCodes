#include <stdio.h>

void compare(TString &dir="./"){
  gROOT->LoadMacro("adamStyles.C");

  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();

  gROOT->LoadMacro("adamCanvasMacros.C");
  gROOT->LoadMacro("adamGetObjMacros.C");
  gROOT->LoadMacro("adamDrawMacros.C");

  TSortedList * directories = new TSortedList();
  directories->Add(getDirectory("TFS_0_500.root","/paramswitch") );
  directories->Add(getDirectory("TFS_1.root","/paramswitch") );
  directories->Add(getDirectory("TFS_5.root","/paramswitch") );
  directories->Add(getDirectory("TFS_10.root","/paramswitch") );
  directories->Add(getDirectory("TFS_100.root","/paramswitch") );
  directories->Add(getDirectory("TFS_200.root","/paramswitch") );
  directories->Add(getDirectory("TFS_500.root","/paramswitch") );
  directories->Add(getDirectory("TFS_1000.root","/paramswitch") );
  directories->Add(getDirectory("TFS_2000.root","/paramswitch") );
  directories->Add(getDirectory("TFS_3000.root","/paramswitch") );

  TString histos[] = {
    //"h_switch_tracker", "h_switch_global"
    "h_switch_50", "h_switch_100", "h_switch_150",  "h_switch_200", "h_switch_250"
    //"h_switch_s1","h_switch_s2","h_switch_s3"
  }


  //  TIter iterH(histos);

  TDirectory* c;
  TString* h;
  TH1 * tmp;
  TH1 * tmp2;
  TH1 * tmp3;
  TH1 * tmp4;
  TH1 * tmp5;
  TH1 * tmp6;
  TH1 * tmp7;
  TGraph * tmpG;
  TGraph * tmpG2;

  int i = 1;


  
  TIter iter(directories);
  while( (c = (TDirectory *)iter()) )
    { 
      //TString iString(i);
      //cout << iString << endl;
      TCanvas * c1_n1 = new TCanvas(); //TString("trk_glb_")+iString,TString("trk_glb_")+iString);
      TLegend * leg = new TLegend(0.75,0.75,0.9,0.9);

      if(histos[0]){
      tmp = getHistogram(c,histos[0]);
      tmp->SetMarkerColor(kRed);
      tmp->SetLineColor(kRed);
      tmp->Draw();
      leg->AddEntry(tmp,histos[0],"l");
      }

      if(histos[1]){
      tmp2 = getHistogram(c,histos[1]);
      tmp2->SetMarkerColor(kBlue);
      tmp2->SetLineColor(kBlue);
      tmp2->Draw("same");
      leg->AddEntry(tmp2,histos[1],"l");
      }
      
      if(histos[2]){
      tmp3 = getHistogram(c,histos[2]);
      tmp3->SetMarkerColor(kBlack);
      tmp3->SetLineColor(kBlack);
      tmp3->Draw("same");
      leg->AddEntry(tmp3,histos[2],"l");
      }
      
      
      if(histos[3]){
      tmp4 = getHistogram(c,histos[3]);
      tmp4->SetMarkerColor(kGreen);
      tmp4->SetLineColor(kGreen);
      tmp4->Draw("same");
      leg->AddEntry(tmp4,histos[3],"l");
      }

      if(histos[4]){
      tmp5 = getHistogram(c,histos[4]);
      tmp5->SetMarkerColor(7);
      tmp5->SetLineColor(7);
      tmp5->Draw("same");
      leg->AddEntry(tmp5,histos[4],"l");
      }
      

      c1_n1->SetLogy(1);
      i++;
      leg->Draw("same");
    }

}
