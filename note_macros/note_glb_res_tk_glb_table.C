#include <stdio.h>

TLegend * note_glb_res_tk_glb_table()
{
  FILE * pFile;
  pFile = fopen ("FigParam/resolutions/global/tab_glb_res_pt_trk_glb_inclusive.tex","w");

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

  TString dirName_("RecoMuon_TrackAssoc/Trk");
  TString dirName_2("RecoMuon_TrackAssoc/Glb");

  TList * fileList = makeFileCollection("my2112FileList_all.txt");

  TString directoriesTrk[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName_,
  }

  TString directoriesGlb[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName_2,
  }

  TList * dirListTrk = makeDirectoryCollection(fileList,directoriesTrk,1);
  TList * dirListGlb = makeDirectoryCollection(fileList,directoriesGlb,1);

  TList * collectionBarrelTrk = makeObjectCollection(dirListTrk,"ErrPt_barrel");
  TList * collectionOverlapTrk = makeObjectCollection(dirListTrk,"ErrPt_overlap");
  TList * collectionEndcapTrk = makeObjectCollection(dirListTrk,"ErrPt_endcap");

  TList * collectionBarrelGlb = makeObjectCollection(dirListGlb,"ErrPt_barrel");
  TList * collectionOverlapGlb = makeObjectCollection(dirListGlb,"ErrPt_overlap");
  TList * collectionEndcapGlb = makeObjectCollection(dirListGlb,"ErrPt_endcap");

  //collectionBarrelTrk->Print();

  double pt[] = {1.,5.,10.,100.,200.,500.,1000.,2000.,3000.};


  TIter iterBTrk(collectionBarrelTrk);
  TH1 * hBTrk;
  TIter iterOTrk(collectionOverlapTrk);
  TH1 * hOTrk;
  TIter iterETrk(collectionEndcapTrk);
  TH1 * hETrk;

  TIter iterBGlb(collectionBarrelGlb);
  TH1 * hBGlb;
  TIter iterOGlb(collectionOverlapGlb);
  TH1 * hOGlb;
  TIter iterEGlb(collectionEndcapGlb);
  TH1 * hEGlb;

  fprintf(pFile,
"\\begin{table}[!hb]
  \\centering
  \\caption[Comparison of the $q/\\pt$ resolution for the tracker system and the combined tracker and muon systems]{Resolution on $q/\\pt$ divided by pseudorapidity regions for
    the different muon reconstruction steps. \\label{tab:qptcomparison}} 
  \\tspace
  \\resizebox*{1\\textwidth}{!}{
    \\begin{tabular}{cccc}
    \\toprule
                        & \\multicolumn{3}{c}{$R(q/\\pt)$ (\\%)}  \\\\
    \\pt{} (\\GeVc)       & \\multicolumn{3}{c}{Tracker/Global} \\\\
                        &  Barrel & Overlap & End-caps \\\\
    \\midrule \n"
	  );

  int i = 0;
  while ( (hBTrk = (TH1 *)iterBTrk() ) )  {
    
    hOTrk = (TH1*)iterOTrk();
    hETrk = (TH1*)iterETrk();

    hBGlb = (TH1*)iterBGlb();
    hOGlb = (TH1*)iterOGlb();
    hEGlb = (TH1*)iterEGlb();

    std::pair<double,double> fitPairBTrk = fit(hBTrk,2);
    std::pair<double,double> fitPairOTrk = fit(hOTrk,2);
    std::pair<double,double> fitPairETrk = fit(hETrk,2);

    std::pair<double,double> fitPairBGlb = fit(hBGlb,2);
    std::pair<double,double> fitPairOGlb = fit(hOGlb,2);
    std::pair<double,double> fitPairEGlb = fit(hEGlb,2);

    fprintf(pFile,"$%d$ & $%7.4f \\pm %4.2e $",pt[i],100*fitPairBTrk.first,100*fitPairBTrk.second);
    fprintf(pFile," / $%7.4f \\pm %4.2e $",100*fitPairBGlb.first,100*fitPairBGlb.second);

    fprintf(pFile," & $%7.4f \\pm %4.2e $",100*fitPairOTrk.first,100*fitPairOTrk.second);
    fprintf(pFile," / $%7.4f \\pm %4.2e $",100*fitPairOGlb.first,100*fitPairOGlb.second);

    fprintf(pFile," & $%7.4f \\pm %4.2e $",100*fitPairETrk.first,100*fitPairETrk.second);
    fprintf(pFile," / $%7.4f \\pm %4.2e $",100*fitPairEGlb.first,100*fitPairEGlb.second);
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
