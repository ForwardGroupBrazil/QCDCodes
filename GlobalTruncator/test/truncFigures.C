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

  TString directories[] = { "/DQMData/RecoMuonV/Trunc/globalTrunc_"} //first1_/"}
  TString name0( "/DQMData/RecoMuonV/Trunc/globalTrunc_");
  TString name("first1_");

  std::vector<TString> newDirectories;
  newDirectories.push_back("first1_/");
  newDirectories.push_back("first2_/");
  newDirectories.push_back("first3_/");
  newDirectories.push_back("first4_/");

  TString directories2[newDirectories.size()];


  cout << "newDirectories size " << newDirectories.size() << endl;

 TList * graphListb = new TList();
 TList * graphListo = new TList();
 TList * graphListe = new TList();

 TString directories3[1];

  int jj = 0;
  for( ; jj < newDirectories.size(); ++jj){
    cout <<"line 40 " << jj <<endl;
    directories3[0] = 0;  cout <<"line 41" <<endl;
    cout << "VECTOR " << newDirectories[jj].Data() << endl;
    cout <<"line34 "<< jj << endl;
    directories3[0] =  TString(name0+newDirectories[jj]);
    cout <<"line36"<<endl;
    TList * dirList2 = makeDirectoryCollection(fileList,directories3,1);

    TIter iter(dirList2);
    TDirectory * c;
    while ( (c= (TDirectory *)iter()) ) {
      cout << "The Directory Collection " << c->GetName() << endl;
    }
    cout <<"line 52" <<endl;
    TList * barrel_collection = makeObjectCollection(dirList2,"ResPtBarrel");
    cout <<"line 54" <<endl;
    /*
    const int nn_b = barrel_collection->GetSize();
    cout <<"line 56 " << nn_b <<endl;
    double *x_b = new double[nn_b];
    double *e_b = new double[nn_b];
    double *yw_b = new double[nn_b];
    double *eyw_b = new double[nn_b];
    cout <<"line 61" <<endl;

  TIter iterGb(barrel_collection);
  TH1* hb;
  int jb =0;    cout <<"line 65 " << jb << endl;
  while( (hb=(TH1*)iterGb()) ) {
    cout <<"line 67 " << jb << endl;
    std::pair<double,double> fitResultb = fit(hb,2);
    yw_b[jb] = fitResultb.first;
    eyw_b[jb] = fitResultb.second;
    x_b[jb] = binCenter[jb];
    e_b[jb] = 0.;
    
    jb++;
  }
    cout <<"line 76" <<endl;
  TGraphErrors *grWb = new TGraphErrors(nn_b,x_b,yw_b,e_b,eyw_b);
  grWb->SetTitle( newDirectories[jj].Data());
  graphListb->Add(grWb);
  cout <<"line 80" <<endl;
    */
    TGraphErrors * grWb = makeTGraphFromTHCollection(barrel_collection,binCenter,2);
    grWb->SetTitle( newDirectories[jj].Data());
    graphListb->Add(grWb);
    cout <<"line 80" <<endl;
  

  TList * overlap_collection = makeObjectCollection(dirList2,"ResPtOverlap");
  /*
  const int nn_o = overlap_collection->GetSize();
  double *x_o = new double[nn_o];
  double *e_o = new double[nn_o];
  double *yw_o = new double[nn_o];
  double *eyw_o = new double[nn_o];
  cout <<"line 90" <<endl;

  TIter iterGo(overlap_collection);
  TH1* ho;
  int jo =0;
  while( (ho=(TH1*)iterGo()) ) {
    cout <<"line 96 " << jo <<endl;
    std::pair<double,double> fitResulto = fit(ho,2);
    yw_o[jo] = fitResulto.first;
    eyw_o[jo] = fitResulto.second;
    x_o[jo] = binCenter[jo];
    e_o[jo] = 0.;
    
    jo++;
  }
  cout <<"line 105" <<endl;
  TGraphErrors *grWo = new TGraphErrors(nn_o,x_o,yw_o,e_o,eyw_o);
  grWo->SetTitle( newDirectories[jj].Data());
  graphListo->Add(grWo);
  cout <<"line 109" <<endl;
  */

    TGraphErrors * grWo = makeTGraphFromTHCollection(overlap_collection,binCenter,2);
    grWo->SetTitle( newDirectories[jj].Data());
    graphListo->Add(grWo);

  TList * endcap_collection = makeObjectCollection(dirList2,"ResPtEndcap");
  /*
  const int nn_e = endcap_collection->GetSize();
  double *x_e = new double[nn_e];
  double *e_e = new double[nn_e];
  double *yw_e = new double[nn_e];
  double *eyw_e = new double[nn_e];
  cout <<"line 116" <<endl;

  TIter iterGe(endcap_collection);
  TH1* he;
  int je =0;
  while( (he=(TH1*)iterGe()) ) {
    cout <<"line 122 " << je <<endl;
    std::pair<double,double> fitResulte = fit(he,2);
    yw_e[je] = fitResulte.first;
    eyw_e[je] = fitResulte.second;
    x_e[je] = binCenter[je];
    e_e[je] = 0.;
    
    je++;
  }
  cout <<"line 131" <<endl;
  TGraphErrors *grWe = new TGraphErrors(nn_e,x_e,yw_e,e_e,eyw_e);
  grWe->SetTitle( newDirectories[jj].Data());
  graphListe->Add(grWe);
  */
    TGraphErrors * grWe = makeTGraphFromTHCollection(endcap_collection,binCenter,2);
    grWe->SetTitle( newDirectories[jj].Data());
    graphListe->Add(grWe);

  cout <<"line 136" <<endl;
    
  }
  new TCanvas();
  drawObjectCollection(graphListb,true);

  new TCanvas();
  drawObjectCollection(graphListo,true);

  new TCanvas();
  drawObjectCollection(graphListe,true);
  
  return;

    cout <<"line38"<<endl;
  TList * dirList = makeDirectoryCollection(fileList,directories,1);

  //  TList * dirList2 = makeDirectoryCollection(fileList,directories2,jj);
    cout <<"line41"<<endl;


  TList * collection = makeObjectCollection(dirList,"ResPtBarrel");

  const int nn = collection->GetSize();
  double *x = new double[nn];
  double *e = new double[nn];
  double *yw = new double[nn];
  double *eyw = new double[nn];


  TIter iterG(collection);
  TH1* h;
  int j =0;
  while( (h=(TH1*)iterG()) ) {
    std::pair<double,double> fitResult = fit(h,2);
    yw[j] = fitResult.first;
    eyw[j] = fitResult.second;
    x[j] = binCenter[j];
    e[j] = 0.;
    
    j++;
  }



  TGraphErrors *grW = new TGraphErrors(nn,x,yw,e,eyw);
  grW->SetTitle(name);

  TList * graphList = new TList();
  graphList->Add(grW);

  drawObjectCollection(graphList,true);


}


TGraph*
makeTGraphFromTHCollection(TList * inlist_, double* binCenter, int type=2 )
{
const int nn_b = inlist_->GetSize();
   double *x_b = new double[nn_b];
    double *e_b = new double[nn_b];
    double *yw_b = new double[nn_b];
    double *eyw_b = new double[nn_b];

  TIter iterGb(inlist_);
  TH1* hb;
  int jb =0;    cout <<"line 65 " << jb << endl;
  while( (hb=(TH1*)iterGb()) ) {
    cout <<"line 67 " << jb << endl;
    std::pair<double,double> fitResultb = fit(hb,2);
    yw_b[jb] = fitResultb.first;
    eyw_b[jb] = fitResultb.second;
    x_b[jb] = binCenter[jb];
    e_b[jb] = 0.;
    
    jb++;
  }

  TGraphErrors *grWb = new TGraphErrors(nn_b,x_b,yw_b,e_b,eyw_b);

  return grWb;

}

