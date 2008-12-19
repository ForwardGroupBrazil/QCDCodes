TLegend * truncFigures()
{
  gROOT->LoadMacro("adamStyles.C");
  
  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();
  
  gROOT->LoadMacro("adamCanvasMacros.C");
  gROOT->LoadMacro("adamDrawMacros.C");
  gROOT->LoadMacro("adamFitMacros.C");
  gROOT->LoadMacro("adamGetObjMacros.C");
  gROOT->LoadMacro("adamMakeCollectionMacros.C");
  
  TList * fileList = makeFileCollection("my220FileList.txt");
  fileList->Print();
  double binCenter[] = {100., 200.,300.,400.};  

  TString name0( "/DQMData/RecoMuonV/Trunc/globalTrunc_");

  std::vector<TString> newDirectories;
  newDirectories.push_back("first1_/");
  newDirectories.push_back("first2_/");
  newDirectories.push_back("first3_/");
  newDirectories.push_back("first4_/");

  cout << "newDirectories size " << newDirectories.size() << endl;

  TList * graphListb = new TList();
  TList * graphListo = new TList();
  TList * graphListe = new TList();

  TString directories[1];
  
  int jj = 0;
  for( ; jj < newDirectories.size(); ++jj){
    directories[0] = 0;
    directories[0] =  TString(name0+newDirectories[jj]);

    TList * dirList = makeDirectoryCollection(fileList,directories,1);

    /*
    TIter iter(dirList);
    TDirectory * c;
    while ( (c= (TDirectory *)iter()) ) {
      cout << "The Directory Collection " << c->GetName() << endl;
    }
    */

    TList * barrel_collection = makeObjectCollection(dirList,"ResPtBarrel");
    TList * overlap_collection = makeObjectCollection(dirList,"ResPtOverlap");
    TList * endcap_collection = makeObjectCollection(dirList,"ResPtEndcap");

    TGraphErrors * grWb = makeTGraphFromTHCollection(barrel_collection,binCenter,2);
    grWb->SetTitle( newDirectories[jj].Data());
    graphListb->Add(grWb);

    TGraphErrors * grWo = makeTGraphFromTHCollection(overlap_collection,binCenter,2);
    grWo->SetTitle( newDirectories[jj].Data());
    graphListo->Add(grWo);

    TGraphErrors * grWe = makeTGraphFromTHCollection(endcap_collection,binCenter,2);
    grWe->SetTitle( newDirectories[jj].Data());
    graphListe->Add(grWe);
  }
  
  newCanvas("c_barrel","Barrel");
  TLegend * b_legend = drawObjectCollection(graphListb,true);
  b_legend->SetHeader("0 < |#eta| < 0.8");
  
  newCanvas("c_overlap","Overlap");
  TLegend * o_legend = drawObjectCollection(graphListo,true);
  o_legend->SetHeader("0.8 < |#eta| < 1.2");
  
  newCanvas("c_endcap","Endcap");
  TLegend * e_legend = drawObjectCollection(graphListe,true);
  e_legend->SetHeader("1.2 < |#eta| < 2.4");
  
  return;
  
  
}


TGraph*
makeTGraphFromTHCollection(TList * inlist_, double* binCenter, int type=2 )
{
  const int nn = inlist_->GetSize();
  double *x = new double[nn];
  double *e = new double[nn];
  double *yw = new double[nn];
  double *eyw = new double[nn];
  
  TIter iterGb(inlist_);
  TH1* hb;
  int jb =0;
  while( (hb=(TH1*)iterGb()) ) {
    std::pair<double,double> fitResultb = fit(hb,2);
    yw[jb] = fitResultb.first;
    eyw[jb] = fitResultb.second;
    x[jb] = binCenter[jb];
    e[jb] = 0.;
    
    jb++;
  }
  
  TGraphErrors *grWb = new TGraphErrors(nn,x,yw,e,eyw);
  
  return grWb;
  
}

