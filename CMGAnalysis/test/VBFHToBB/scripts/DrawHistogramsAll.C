#include "DrawHistograms.C"
void DrawHistogramsAll()
{
  const int N = 34;
  TString HISTO[N]  = {"hMet","hMbbCor","hMqq","hdPhibb","hdEtaqq","hBDT","hMLP","hMbbCorCut","hJetEtaBtag0","hJetEtaBtag1","hJetEtaBtag2","hJetEtaBtag3",
"hJetPt0","hJetPt1","hJetPt2","hJetPt3","hJetPtBtag0","hJetPtBtag1","hJetPtBtag2","hJetPtBtag3","hJetBtag0","hJetBtag1",
"hJetQGL2","hJetQGL3","hRho","hdEtaqqDiff","hEtaBoost","hsoftHt","hsoftMulti","hJetPt4Mqq","hJetPuMvaBtag0","hJetPuMvaBtag1",
"hJetPuMvaBtag2","hJetPuMvaBtag3"};
  TString XTITLE[N] = {"MET (GeV)","Regressed M_{bb} (GeV)","M_{qq} (GeV)","#Delta#phi_{bb}","|#Delta#eta_{qq}|","BDT Output","ANN Output","Regressed M_{bb} (GeV)","Jet0 #eta","Jet1 #eta","Jet2 #eta","Jet3 #eta","Jet0 p_{T} (GeV)","Jet1 p_{T} (GeV)","Jet2 p_{T} (GeV)","Jet3 p_{T} (GeV)","BJet0 p_{T} (GeV)","BJet1 p_{T} (GeV)","BJet2 p_{T} (GeV)","BJet3 p_{T} (GeV)","btag0","btag1","qgl2","qgl3","rho","dEtaMax - dEtaqq","EtaBoostqq","Soft H_{T} (GeV)","Soft Multi","Jet P_{T,5}/mqq",
"puMvaBtag0","puMvaBtag1","puMvaBtag2","puMvaBtag3"};
  float YMAX[N] = {1e+6,1e+6,1e+6,1e+6,1e+6,1e+6,1e+6,1e+5,1e+6,1e+6,1e+6,1e+6,1e+6,1e+6,1e+6,1e+6,1e+6,1e+6,1e+6,1e+6,2e+6,1e+6,1e+6,1e+6,1e+6,1e+7,1e+6,1e+6,1e+6,1e+6,1e+7,1e+7,1e+7,1e+7};
  float XMIN[N] = {0,0,300,0,2.5,-0.4,-0.4,0,-5,-5,-5,-5,85,70,60,40,40,40,40,40,,0,0,0,0,0,0,-4,0,0,0,-1,-1,-1,-1};
  float XMAX[N] = {500,3000,4000,3.2,9,0.4,1.4,500,5,5,5,5,1000,800,500,350,1000,800,800,500,1,1,1,1,40,3,4,200,50,1,1,1,1,1};
  bool LOGY[N]  = {true,true,true,false,true,true,true,true,false,false,false,false,true,true,true,true,true,true,true,true,true,false,false,false,true,
true,false,false,false,false,false,false,false,false};
  int REBIN[N] = {1,2,2,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,2,2,2,2,2};
  for(int i=0;i<34;i++) {
    DrawHistograms(HISTO[i],XTITLE[i],YMAX[i],XMIN[i],XMAX[i],REBIN[i],true,true);
  }
}
