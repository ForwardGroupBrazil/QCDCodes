#include <stdio.h>

void note_tk_eff_table()
{
  FILE * pFile;
  //  pFile = fopen ("FigGLB/efficiencies/tracker/tab_tk_eff_inclusive.tex","w");
  pFile = fopen ("FigGLB/efficiencies/trackerMuons/tab_tkmu_eff_inclusive.tex","w");

  gROOT->LoadMacro("adamStyles.C");
  
  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();
  
  //  gROOT->LoadMacro("adamCanvasMacros.C");
  //  gROOT->LoadMacro("adamDrawMacros.C");
  //  gROOT->LoadMacro("adamFitMacros.C");
  gROOT->LoadMacro("adamGetObjMacros.C");
  gROOT->LoadMacro("adamMakeCollectionMacros.C");
  
  TList * fileList = makeFileCollection("my2112FileList_high.txt");
  
  fileList->Print();

  TString directories[] = {
    //    "/DQMData/Run 1/RecoMuonV/Run summary/MultiTrack/general_tpToTkmuAssociation"
    "/DQMData/Run 1/RecoMuonV/Run summary/MultiTrack/good_TrackerOnly_tpToTkmutkAssociation"
  }
  TList * dirList = makeDirectoryCollection(fileList,directories,1);

  TIter iter(dirList);
  TDirectory * c;
  while ( (c= (TDirectory *)iter()) ) {
    cout << "The Directory Collection " << c->GetName() << endl;
  }

  TList * collection = makeObjectCollection(dirList,"efficPt");
  collection->Print();

  double pt[] = {5.,10.,100.,200.,500.,1000.,2000.,3000.};

  TIter iterG(collection);
  TH1* h;
  int i=0;
  fprintf(pFile,"
\\begin{table}[!h]
  \\centering
  \\caption{Integrated efficiencies for the tracker tracks
  reconstruction. \\label{tab:tk_eff}}
  \\tspace
  \\begin{tabular}{cccc}
    \\toprule
    \\pt{} sample             & $\\varepsilon_{tk}$ 
    \\\\
    (\\GeVc)                  & (\\%)  \\\\

    \\midrule

");

  while( (h=(TH1*)iterG()) ) {

    int bin = h->Fill(pt[i],0);
    fprintf(pFile,"$%d$ & $%7.4f \\pm %4.2e $",pt[i],100*h->GetBinContent(bin),100*h->GetBinError(bin));
    fprintf(pFile,"\\\\ \n");
    i++;
  }

  fprintf(pFile,"    \\bottomrule
  \\end{tabular}
\\end{table}");
  
  fclose(pFile);

}
