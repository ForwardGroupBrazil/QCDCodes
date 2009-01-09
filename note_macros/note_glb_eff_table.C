#include <stdio.h>

void note_glb_eff_table()
{
  FILE * pFile;
  pFile = fopen ("tab_glb_eff_inclusive.tex","w");

  gROOT->LoadMacro("adamStyles.C");
  
  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();
  
  //  gROOT->LoadMacro("adamCanvasMacros.C");
  //  gROOT->LoadMacro("adamDrawMacros.C");
  //  gROOT->LoadMacro("adamFitMacros.C");
  gROOT->LoadMacro("adamGetObjMacros.C");
  gROOT->LoadMacro("adamMakeCollectionMacros.C");
  
  TList * fileList = makeFileCollection("my2112FileList2.txt");
  
  fileList->Print();

  TString directories[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/MultiTrack/globalMuons_tpToGlbAssociation"
  }
  TList * dirList = makeDirectoryCollection(fileList,directories,1);

  TString directoriesComp[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/MultiTrack"
  }
  TList * dirListComp = makeDirectoryCollection(fileList,directoriesComp,1);

  TIter iter(dirList);
  TDirectory * c;
  while ( (c= (TDirectory *)iter()) ) {
    cout << "The Directory Collection " << c->GetName() << endl;
  }

  //  TList * collection = makeObjectCollection(dirList,"ptres_vs_eta_Sigma");
  TList * collectionCompSta = makeObjectCollection(dirListComp,"Eff_GlbSta_Pt");
  collectionCompSta->Print();

  TList * collectionCompTk = makeObjectCollection(dirListComp,"Eff_GlbTk_Pt");
  collectionCompTk->Print();

  TList * collection = makeObjectCollection(dirList,"efficPt");
  collection->Print();

  double pt[] = {10.,100.,1000.};

  TIter iterG(collection);
  TH1* h;
  TIter iterCompSta(collectionCompSta);
  TH1* hSta;
  TIter iterCompTk(collectionCompTk);
  TH1* hTk;
  int i=0;
  fprintf(pFile,"\\begin{table}[!h]
  \\centering
  \\caption{Integrated efficiencies for the global muon
  reconstruction. \\label{tab:glb_eff}}
  \\tspace
  \\begin{tabular}{ccccc}
    \\toprule
    \\pt{} sample         &  $\\varepsilon_{glb}$  
                         &  $\\varepsilon_{glb,sa}$ 
                         &  $\\varepsilon_{glb,tk}$ 
                          \\\\ 
    (\\GeVc) & (\\%) & (\\%) & (\\%)  \\\\
	  \\midrule
");
  while( (h=(TH1*)iterG()) ) {

    hSta = (TH1*)iterCompSta();
    hTk = (TH1*)iterCompTk();

    int bin = h->Fill(pt[i],0);
    int binSta = hSta->Fill(pt[i],0);
    int binTk = hTk->Fill(pt[i],0);
    fprintf(pFile,"$%d$ & $%7.4f \\pm %4.2e $",pt[i],100*h->GetBinContent(bin),100*h->GetBinError(bin));
    fprintf(pFile," & $ %7.4f \\pm %4.2e $",100*hSta->GetBinContent(binSta),100*hSta->GetBinError(binSta));
    fprintf(pFile," & $ %7.4f \\pm %4.2e $",100*hTk->GetBinContent(binTk),100*hTk->GetBinError(binTk));
    fprintf(pFile,"\\\\ \n");
    i++;
  }

  fprintf(pFile,"    \\bottomrule
  \\end{tabular}
\\end{table}");
  
  fclose(pFile);

}
