#include <stdio.h>

TLegend * note_sta_res_noUp_table()
{
  FILE * pFile;
  //pFile = fopen ("tmpStaTable.tex","w");
  pFile = fopen ("FigParam/resolutions/standAlone/tab_sta_res_ptb_up_noUp_inclusive.tex","w");
  fprintf(pFile,
"\\begin{table}[!hb]
  \\centering
   \\caption[Stand-alone muon $\\pt$ resolution before and after the constraint at the vertex]{Stand-alone muon resolution on $\\pt$ before and after the constraint at
     the vertex, divided by pseudorapidity regions. \\label{tab:ptb_constraintAtVtx}} 
   \\tspace
   \\resizebox*{1\\textwidth}{!}{
     \\begin{tabular}{cccc}
       \\toprule
       & \\multicolumn{3}{c}{$R(q/\\pt)$ (\\%)}  \\\\
       \\pt{} (\\GeVc)       & \\multicolumn{3}{c}{Before/after the vertex constraint} \\\\
       &  Barrel & Overlap & End-caps \\\\
       \\midrule \n"
	  );
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

  TList * fileList = makeFileCollection("my2112FileList_all.txt");

  double pt[] = {1.,5.,10.,100.,200.,500.,1000.,2000.,3000.};

  TString dirName_("MultiTrack/standAloneMuons_UpdatedAtVtx_tpToStaAssociation");
  TString dirName2_("MultiTrack/standAloneMuons_tpToStaNoUpAssociation");

  TString dirName3_("RecoMuon_TrackAssoc/Sta");
  TString dirName4_("RecoMuon_NoUpTrackAssoc/Sta");

  TString directories[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName_
  }

  TString directories2[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName2_
  }

  TString directoriesRM3[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName3_
  }

  TString directoriesRM4[] = {
    "/DQMData/Run 1/RecoMuonV/Run summary/"+dirName4_
  }

  TList * dirList = makeDirectoryCollection(fileList,directories,1);
  TList * dirList2 = makeDirectoryCollection(fileList,directories2,1);
  TList * dirList3 = makeDirectoryCollection(fileList,directoriesRM3,1);
  TList * dirList4 = makeDirectoryCollection(fileList,directoriesRM4,1);

  //----------------------------------------------------
  // Now do each specific figure
  //----------------------------------------------------

  //----------------------------------------------------
  // Resolutions
  //----------------------------------------------------
  
  TString histo[] = {
		     "ptres",
		     "ErrP",
		     "ErrQPt",
		     "ErrPt"};

  int jj = 2;

  TList * glbCol;
  TList * trkCol;

  if(jj==0){
    glbCol = makeObjectCollection(dirList,"ptres_vs_eta");
    trkCol = makeObjectCollection(dirList2,"ptres_vs_eta");
  }

  if(jj==1) {
    glbCol = makeObjectCollection(dirList3,"ErrP_vs_Eta");
    trkCol = makeObjectCollection(dirList4,"ErrP_vs_Eta");
  }

  if(jj==2) {
    glbCol = makeObjectCollection(dirList3,"ErrQPt_vs_Eta");
    trkCol = makeObjectCollection(dirList4,"ErrQPt_vs_Eta");
  }

  if(jj==3) {
    glbCol = makeObjectCollection(dirList3,"ErrPt_vs_Eta");
    trkCol = makeObjectCollection(dirList4,"ErrPt_vs_Eta");
  }

  TList * fitColGlb = makeFitCollection(glbCol,2);
  TList * fitColTrk = makeFitCollection(trkCol,2);

  TList * trackerBuffer_B = makeFitCollection(fitColTrk,1);
  TList * trackerBuffer_O = makeFitCollection(fitColTrk,2);
  TList * trackerBuffer_E = makeFitCollection(fitColTrk,3);

  TList * globalBuffer_B = makeFitCollection(fitColGlb,1);
  TList * globalBuffer_O = makeFitCollection(fitColGlb,2);
  TList * globalBuffer_E = makeFitCollection(fitColGlb,3);

  trackerBuffer_B->Print();
  globalBuffer_B->Print();
  
  TIter tkIter_B(trackerBuffer_B);
  TIter glbIter_B(globalBuffer_B);
  TIter tkIter_O(trackerBuffer_O);
  TIter glbIter_O(globalBuffer_O);
  TIter tkIter_E(trackerBuffer_E);
  TIter glbIter_E(globalBuffer_E);

  TObjString * Tk_B;
  TObjString * Glb_B;
  TObjString * Tk_O;
  TObjString * Glb_O;
  TObjString * Tk_E;
  TObjString * Glb_E;

  int i = 0;

  while( (Tk_B=(TObjString*)tkIter_B())) {
    Glb_B = (TObjString*)glbIter_B();
    Tk_O = (TObjString*)tkIter_O();
    Glb_O = (TObjString*)glbIter_O();
    Tk_E = (TObjString*)tkIter_E();
    Glb_E = (TObjString*)glbIter_E();
    
    fprintf(pFile,"$%d$ & %s / %s",pt[i],Tk_B->GetString().Data(),Glb_B->GetString().Data());
    fprintf(pFile," & %s / %s",Tk_O->GetString().Data(),Glb_O->GetString().Data());
    fprintf(pFile," & %s / %s",Tk_E->GetString().Data(),Glb_E->GetString().Data());

   fprintf(pFile,"\\\\ \n");
    
    i++;
    
  }

  fprintf(pFile,"    \\bottomrule
  \\end{tabular}
}
\\end{table}");
  
  fclose(pFile);  
}
