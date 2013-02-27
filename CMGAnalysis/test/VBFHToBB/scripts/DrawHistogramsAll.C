#include "DrawHistograms.C"
void DrawHistogramsAll()
{
  const int N = 35;

  TString HISTO[N]  = {"hnVtx","hcosTheta","hcosAlpha","hmet","hmbbCor","hmqq","hdPhibb","hdEtaqq","hMLP","hmbbCorCut","hjetEtaBtag0","hjetEtaBtag1","hjetEtaBtag2","hjetEtaBtag3",
"hjetPt0","hjetPt1","hjetPt2","hjetPt3","hjetPtBtag0","hjetPtBtag1","hjetPtBtag2","hjetPtBtag3","hjetBtag0","hjetBtag1",
"hjetQGL2","hjetQGL3","hrho","hdEtaqqDiff","hetaBoostqq","hsoftHt","hsoftMulti","hjetPuMva0","hjetPuMva1",
"hjetPuMva2","hjetPuMva3"};
  TString XTITLE[N] = {"Number of vertices","|cos#theta|","|cos#alpha|","MET (GeV)","Regressed M_{bb} (GeV)","M_{qq} (GeV)","#Delta#phi_{bb}","|#Delta#eta_{qq}|","ANN Output","Regressed M_{bb} (GeV)","Jet0 #eta","Jet1 #eta","Jet2 #eta","Jet3 #eta","Jet0 p_{T} (GeV)","Jet1 p_{T} (GeV)","Jet2 p_{T} (GeV)","Jet3 p_{T} (GeV)","BJet0 p_{T} (GeV)","BJet1 p_{T} (GeV)","BJet2 p_{T} (GeV)","BJet3 p_{T} (GeV)","btag0","btag1","qgl2","qgl3","rho","dEtaMax - dEtaqq","EtaBoostqq","Soft H_{T} (GeV)","Soft Multi",
"puMvaBtag0","puMvaBtag1","puMvaBtag2","puMvaBtag3"};
  float YMAX[N] = {1e+6,1e+6,1e+6,1e+6,1e+6,1e+6,1e+6,1e+6,1e+6,1e+5,1e+6,1e+6,1e+6,1e+7,1e+7,1e+7,1e+7,1e+7,1e+6,1e+6,1e+6,1e+6,1e+7,1e+7,1e+7,1e+7,1e+6,1e+7,1e+6,1e+6,1e+6,1e+7,1e+7,1e+7,1e+7};
  float XMIN[N] = {0,-1,-1,0,0,300,0,2.5,-0.4,0,-5,-5,-5,-5,85,70,60,40,40,40,40,40,0,0,0,0,0,0,-4,0,0,-1,-1,-1,-1};
  float XMAX[N] = {50,1,1,500,3000,4000,1.99,9,1.4,500,5,5,5,5,1000,800,500,350,1000,800,800,500,1,1,1,1,40,3,4,200,50,1,1,1,1};
  bool LOGY[N]  = {true,false,false,true,true,true,false,true,true,true,false,false,false,false,true,true,true,true,true,true,true,true,true,false,false,false,true,
true,false,false,false,false,false,false,false};
  int REBIN[N] = {1,1,1,1,2,2,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,2,2,2,2};
  for(int i=0;i<N;i++) {
    DrawHistograms(HISTO[i],XTITLE[i],YMAX[i],XMIN[i],XMAX[i],REBIN[i],true,true);
  }
}
