#include <stdio.h>

TLegend * note_sta_res_noUp_table()
{
  FILE * pFile;
  pFile = fopen ("FigParam/resolutions/standAlone/tab_sta_res_pt_up_noUp_inclusive.tex","w");

  //----------------------------------------------------
  // General macro setup
  //----------------------------------------------------

  gROOT->LoadMacro("adamStyles.C");
  
  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();
  
  gROOT->LoadMacro("adamCanvasMacros.C");
  gROOT->LoadMacro("adamDrawMacros.C");
  gROOT->LoadMacro("adamFitMacros.C");
  gROOT->LoadMacro("adamGetObjMacros.C");
  gROOT->LoadMacro("adamMakeCollectionMacros.C");
  
  //----------------------------------------------------
  // General setup for files and directories
  //----------------------------------------------------

  TString dirName_("RecoMuon_NoUpTrackAssoc/Sta");
  TString dirName_2("RecoMuon_TrackAssoc/Sta");
  TString figDirName_("standAlone/staNoUp");
  TString figLongName_("Stand Alone Muons");

  TList * fileList = makeFileCollection("my2112FileList_all.txt");

  TString directoriesNoUp[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName_,
  }

  TString directoriesUp[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName_2,
  }

  TList * dirListNoUp = makeDirectoryCollection(fileList,directoriesNoUp,1);
  TList * dirListUp = makeDirectoryCollection(fileList,directoriesUp,1);

  TList * collectionBarrelNoUp = makeObjectCollection(dirListNoUp,"ErrPt_barrel");
  TList * collectionOverlapNoUp = makeObjectCollection(dirListNoUp,"ErrPt_overlap");
  TList * collectionEndcapNoUp = makeObjectCollection(dirListNoUp,"ErrPt_endcap");

  TList * collectionBarrelUp = makeObjectCollection(dirListUp,"ErrPt_barrel");
  TList * collectionOverlapUp = makeObjectCollection(dirListUp,"ErrPt_overlap");
  TList * collectionEndcapUp = makeObjectCollection(dirListUp,"ErrPt_endcap");

  //collectionBarrelNoUp->Print();

  double pt[] = {1.,5.,10.,100.,200.,500.,1000.,2000.,3000.};


  TIter iterBNoUp(collectionBarrelNoUp);
  TH1 * hBNoUp;
  TIter iterONoUp(collectionOverlapNoUp);
  TH1 * hONoUp;
  TIter iterENoUp(collectionEndcapNoUp);
  TH1 * hENoUp;

  TIter iterBUp(collectionBarrelUp);
  TH1 * hBUp;
  TIter iterOUp(collectionOverlapUp);
  TH1 * hOUp;
  TIter iterEUp(collectionEndcapUp);
  TH1 * hEUp;

  fprintf(pFile,"\\begin{table}[thb]
   \\centering
   \\caption[Stand-alone muon $q/\\pt$ resolution before and after the constraint at vertex]{stand-alone muon resolution on $q/\\pt$ before and after the constraint at
     vertex, divided by pseudorapidity regions. \\label{tab:constraintAtVtx}} 
   \\tspace
   \\resizebox*{1\\textwidth}{!}{
     \\begin{tabular}{cccc}
       \\toprule
       & \\multicolumn{3}{c}{$R(q/\\pt)$ (\\%)}  \\\\
       \\pt{} (\\GeVc)       & \\multicolumn{3}{c}{Before/after the vertex constraint} \\\\
       &  Barrel & Overlap & End-caps \\\\
       \\midrule \n");

  int i = 0;
  while ( (hBNoUp = (TH1 *)iterBNoUp() ) )  {
    
    hONoUp = (TH1*)iterONoUp();
    hENoUp = (TH1*)iterENoUp();

    hBUp = (TH1*)iterBUp();
    hOUp = (TH1*)iterOUp();
    hEUp = (TH1*)iterEUp();

    std::pair<double,double> fitPairBNoUp = fit(hBNoUp,2);
    std::pair<double,double> fitPairONoUp = fit(hONoUp,2);
    std::pair<double,double> fitPairENoUp = fit(hENoUp,2);

    std::pair<double,double> fitPairBUp = fit(hBUp,2);
    std::pair<double,double> fitPairOUp = fit(hOUp,2);
    std::pair<double,double> fitPairEUp = fit(hEUp,2);

    fprintf(pFile,"$%d$ & $%7.4f \\pm %4.2e $",pt[i],100*fitPairBNoUp.first,100*fitPairBNoUp.second);
    fprintf(pFile," / $%7.4f \\pm %4.2e $",100*fitPairBUp.first,100*fitPairBUp.second);

    fprintf(pFile," & $%7.4f \\pm %4.2e $",100*fitPairONoUp.first,100*fitPairONoUp.second);
    fprintf(pFile," / $%7.4f \\pm %4.2e $",100*fitPairOUp.first,100*fitPairOUp.second);

    fprintf(pFile," & $%7.4f \\pm %4.2e $",100*fitPairENoUp.first,100*fitPairENoUp.second);
    fprintf(pFile," / $%7.4f \\pm %4.2e $",100*fitPairEUp.first,100*fitPairEUp.second);
    fprintf(pFile,"\\\\ \n");
    
    i++;
  }

  fprintf(pFile,"    \\bottomrule
  \\end{tabular}
}
\\end{table}");

    fclose(pFile);
  //printCanvasesType(".eps");
}
