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
  double binCenter[] = {10., 100., 200., 500., 1000., 2000., 3000.};  

  TString name0( "/DQMData/RecoMuonV/Trunc/globalTrunc_");

  std::vector<TString> newDirectories;
  /*
  newDirectories.push_back("first1_/");
  newDirectories.push_back("first2_/");
  newDirectories.push_back("first3_/");
  newDirectories.push_back("first4_/");
  */
  /*
  newDirectories.push_back("station1_/");
  newDirectories.push_back("station2_/");
  newDirectories.push_back("station3_/");
  newDirectories.push_back("station4_/");
  */
    
  newDirectories.push_back("picky1_/");
  newDirectories.push_back("picky2_/");
  newDirectories.push_back("picky3_/");
  newDirectories.push_back("picky4_/");
  

  TString labels[] = {"Tracker + 1","Tracker + 2","Tracker + 3","Tracker + 4"}
  //TString labels[] = {"Tracker + Mu1","Tracker + Mu2","Tracker + Mu3","Tracker + Mu4"}

  cout << "newDirectories size " << newDirectories.size() << endl;

  TList * graphListb = new TList();
  TList * graphListo = new TList();
  TList * graphListe = new TList();

  TList * graphListbp = new TList();
  TList * graphListop = new TList();
  TList * graphListep = new TList();

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
    barrel_collection->SetName("barrel"+newDirectories[jj].Chop().Chop());

    TList * overlap_collection = makeObjectCollection(dirList,"ResPtOverlap");
    overlap_collection->SetName("overlap"+newDirectories[jj]);

    TList * endcap_collection = makeObjectCollection(dirList,"ResPtEndcap");
    endcap_collection->SetName("endcap"+newDirectories[jj]);

    TList * barrel_p_collection = makeObjectCollection(dirList,"ResPBarrel");
    barrel_p_collection->SetName("barrel_p_"+newDirectories[jj]);
    TList * overlap_p_collection = makeObjectCollection(dirList,"ResPOverlap");
    overlap_p_collection->SetName("overlap_p_"+newDirectories[jj]);
    TList * endcap_p_collection = makeObjectCollection(dirList,"ResPEndcap");
    endcap_p_collection->SetName("endcap_p_"+newDirectories[jj]);

    TGraphErrors * grWb = makeTGraphFromTHCollection(barrel_collection,binCenter,2);
    grWb->SetTitle( newDirectories[jj].Data());
    graphListb->Add(grWb);

    TGraphErrors * grWo = makeTGraphFromTHCollection(overlap_collection,binCenter,2);
    grWo->SetTitle( newDirectories[jj].Data());
    graphListo->Add(grWo);

    TGraphErrors * grWe = makeTGraphFromTHCollection(endcap_collection,binCenter,2);
    grWe->SetTitle( newDirectories[jj].Data());
    graphListe->Add(grWe);
    /////
    TGraphErrors * grWbp = makeTGraphFromTHCollection(barrel_p_collection,binCenter,2);
    grWbp->SetTitle( newDirectories[jj].Data());
    graphListbp->Add(grWbp);

    TGraphErrors * grWop = makeTGraphFromTHCollection(overlap_p_collection,binCenter,2);
    grWop->SetTitle( newDirectories[jj].Data());
    graphListop->Add(grWop);

    TGraphErrors * grWep = makeTGraphFromTHCollection(endcap_p_collection,binCenter,2);
    grWep->SetTitle( newDirectories[jj].Data());
    graphListep->Add(grWep);
  }
  
  TCanvas * c_barrel = newCanvas("c_barrel","Barrel");
  TLegend * b_legend = drawObjectCollection(graphListb,true,labels);
  b_legend->SetHeader("0 < |#eta| < 0.8");
  c_barrel->SetLogy(1);
  c_barrel->SetLogx(1);
  ((TGraph*)graphListb->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.01,1.);
  
  TCanvas * c_overlap = newCanvas("c_overlap","Overlap");
  TLegend * o_legend = drawObjectCollection(graphListo,true,labels);
  o_legend->SetHeader("0.8 < |#eta| < 1.2");
  c_overlap->SetLogy(1);
  c_overlap->SetLogx(1);
  ((TGraph*)graphListo->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.01,1.);

  TCanvas * c_endcap = newCanvas("c_endcap","Endcap");
  TLegend * e_legend = drawObjectCollection(graphListe,true,labels);
  e_legend->SetHeader("1.2 < |#eta| < 2.4");
  c_endcap->SetLogy(1);
  c_endcap->SetLogx(1);
  ((TGraph*)graphListe->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.01,1.);

  //////

  TCanvas * c_barrel_p = newCanvas("c_barrel_p","Barrel_p");
  TLegend * b_legend_p = drawObjectCollection(graphListbp,true,labels);
  b_legend_p->SetHeader("0 < |#eta| < 0.8");
  c_barrel_p->SetLogy(1);
  c_barrel_p->SetLogx(1);
  ((TGraph*)graphListbp->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.01,1.);
  
  TCanvas * c_overlap_p = newCanvas("c_overlap_p","Overlap_p");
  TLegend * o_legend_p = drawObjectCollection(graphListop,true,labels);
  o_legend_p->SetHeader("0.8 < |#eta| < 1.2");
  c_overlap_p->SetLogy(1);
  c_overlap_p->SetLogx(1);
  ((TGraph*)graphListop->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.01,1.);

  TCanvas * c_endcap_p = newCanvas("c_endcap_p","Endcap_p");
  TLegend * e_legend_p = drawObjectCollection(graphListep,true,labels);
  e_legend_p->SetHeader("1.2 < |#eta| < 2.4");
  c_endcap_p->SetLogy(1);
  c_endcap_p->SetLogx(1);
  ((TGraph*)graphListep->First())->GetHistogram()->GetYaxis()->SetRangeUser(0.01,1.);


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
  
  TCanvas * fitCanvas = newCanvas(inlist_->GetName(),inlist_->GetName(),3,3);
  cout << "TName " << inlist_->GetName() << endl;

  TIter iterGb(inlist_);
  TH1* hb;
  int jb =0;
  int pad = 1;
  while( (hb=(TH1*)iterGb()) ) {
    fitCanvas->cd(pad);
    hb->Draw();
    std::pair<double,double> fitResultb = fit(hb,2);
    yw[jb] = fitResultb.first;
    eyw[jb] = fitResultb.second;
    x[jb] = binCenter[jb];
    e[jb] = 0.;
    //cout << jb << " x " << x[jb] <<" yw " <<  yw[jb] << endl; 
   

    fitCanvas->cd(jb+1);
    hb->Draw();
    
    jb++;
    pad++;
  }
  
  TGraphErrors *grWb = new TGraphErrors(nn,x,yw,e,eyw);
  
  return grWb;
  
}

