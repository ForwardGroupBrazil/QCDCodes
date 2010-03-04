TGraphErrors *
convertTH1toTGraph(TH1* histo_)
{
  TGraphErrors * myGraph = new TGraphErrors(histo_);
  myGraph->SetLineColor(histo_->GetLineColor());
  myGraph->SetLineStyle(histo_->GetLineStyle());
  myGraph->SetMarkerStyle(histo_->GetMarkerStyle());
  return myGraph;
}

TH1* distributionToEff(TH1* dis, bool greaterThanCut) {

  std::ostringstream pprint;
  pprint.str("");

  pprint<<dis->GetName()<<"_integEff";

  TH1* _h2 = dis->Clone(pprint.str().c_str());
  _h2->Sumw2();

   if ( greaterThanCut ) {
      for (int ii = 0; ii != dis->GetNbinsX()+1; ii++) {
        _h2->SetBinContent(ii, dis->Integral(ii, dis->GetNbinsX()+1));
      }
      float intge = dis->Integral(0, dis->GetNbinsX()+1);
      if ( intge > 1e-6)  _h2->Scale(1.0/intge );

   } else {
      for (int ii = 0; ii != dis->GetNbinsX()+1; ii++) {
        _h2->SetBinContent(ii, dis->Integral(0, ii));
      }
      float intge = dis->Integral(0, dis->GetNbinsX());
      if ( intge > 1e-6)  _h2->Scale(1.0/intge );

   }
   return _h2;

}


TLegend *
drawObjectCollection(TList * objectList_, bool draw_=true, TString * legend_ = 0)
{

  //int color[] = { 1, 2, 3, 4, 6, 8, 9,11,12,13,14,15,16,17,18,19,20,21};
  int color[] = { 1, 2, 4, 904, 419, 9,11,12,13,14,15,16,17,18,19,20,21};
  int style[] = {21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38};
  int lstyle[] ={ 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  int type = 0;

  TLegend * myLegend = new TLegend(.35, .15, 0.65, .35, "");
  myLegend->SetFillColor(0);

  TIter hIter(objectList_);
  TObject * obj;
  int i = -1;
  float numerator = 1.;
  while ( (obj=(TObject*)hIter()) ) {
    i++;
    if( obj->IsA()->InheritsFrom("TH1") ) type = 1;
    if( obj->IsA()->InheritsFrom("TH2") ) type = 2;
    if( obj->IsA()->InheritsFrom("TGraph") ) type = 3;

    if(i==0 && type == 1) {
      TH1 * h1 = (TH1*)obj;
      numerator = h1->Integral();
    }

    TString legEntry = (legend_) ? legend_[i] : obj->GetTitle();

    switch(type){
    case 1: {
      TString opt1(color[i]==1 ? "" : "same");
      //TString opt2 = ("E1X0") + opt1;
      TString opt2 = (" PL ") + opt1;
      TH1 * h1 = (TH1*)obj;
      float denominator = (h1->Integral() > 0) ? h1->Integral() : 1.;
      //aaa h1->Scale(numerator/denominator);
      h1->SetLineColor(color[i]);
      h1->SetMarkerColor(color[i]);
      h1->SetMarkerStyle(style[i]);
      h1->SetLineStyle(lstyle[i]);
      if(draw_) h1->Draw(opt2);
      if(draw_) myLegend->AddEntry(h1,legEntry,"PL");
      break;
    }
    case 2: {
      TString opt1(color[i]==1 ? "" : "same");
      TString opt2 = ("E1X0") + opt1;
      TH2 * h2 = (TH2*)obj;
      h2->SetLineColor(color[i]);
      h2->SetMarkerColor(color[i]);
      h2->SetMarkerStyle(style[i]);
      h2->SetLineStyle(lstyle[i]);
      if(draw_) h2->Draw(opt2);
      if(draw_) myLegend->AddEntry(h2,legEntry,"PL");
      break;
    } 
    case 3: {
      TString opt1(color[i]==1 ? "APL" : "PL same");
      TString opt2 = (" ") + opt1;
      TGraph * h3 = (TGraph*)obj;
      h3->SetLineColor(color[i]);
      h3->SetMarkerColor(color[i]);
      h3->SetMarkerStyle(style[i]);
      h3->SetLineStyle(lstyle[i]);
      if(draw_) h3->Draw(opt2);
      if(draw_) myLegend->AddEntry(h3,legEntry,"PL");
      break;
    }
    default: {
      TString opt1(color[i]==1 ? "" : "same");
      TString opt2 = ("E1X0") + opt1;
      TH1 * h4 = (TH1*)obj;
      h4->SetLineColor(color[i]);
      h4->SetMarkerColor(color[i]);
      h4->SetMarkerStyle(style[i]);
      h4->SetLineStyle(lstyle[i]);
      if(draw_) h4->Draw(opt2);
      if(draw_) myLegend->AddEntry(h4,legEntry,"PL");
      break;
    }
    }
  }
  //  if(draw_) myLegend->Draw("same");
  
  return myLegend;
}
