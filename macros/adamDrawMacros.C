TGraphErrors *
convertTH1toTGraph(TH1* histo_)
{
  TGraphErrors * myGraph = new TGraphErrors(histo_);
  myGraph->SetLineColor(histo_->GetLineColor());
  myGraph->SetLineStyle(histo_->GetLineStyle());
  myGraph->SetMarkerStyle(histo_->GetMarkerStyle());
  return myGraph;
}

void
drawObjectCollection(TSortedList * objectList_, bool draw_=true)
{

  int color[] = {1,2,3,4,5,6,7,8,9,10};
  int style[] = {22,23,24,25,26,27,22,23,24,25}
  int type = 0;

  TIter hIter(objectList_);
  TObject * obj;
  int i = -1;
  while ( (obj=(TObject*)hIter()) ) {
    i++;
    if( obj->IsA()->InheritsFrom("TH1") ) type = 1;
    if( obj->IsA()->InheritsFrom("TH2") ) type = 2;
    if( obj->IsA()->InheritsFrom("TGraph") ) type = 3;

    switch(type){
    case 1: {
      TString opt1(color[i]==1 ? "" : "same");
      TString opt2 = ("E1X0") + opt1;
      TH1 * h1 = (TH1*)obj;
      h1->SetLineColor(color[i]);
      h1->SetMarkerColor(color[i]);
      h1->SetMarkerStyle(style[i]);
      if(draw_) h1->Draw(opt2);
      break;
    }
    case 2: {
      TString opt1(color[i]==1 ? "" : "same");
      TString opt2 = ("E1X0") + opt1;
      TH2 * h2 = (TH2*)obj;
      h2->SetLineColor(color[i]);
      h2->SetMarkerColor(color[i]);
      h2->SetMarkerStyle(style[i]);
      if(draw_) h2->Draw(opt2);
      break;
    } 
    case 3: {
      TString opt1(color[i]==1 ? "APL" : "PL same");
      TString opt2 = (" ") + opt1;
      TGraph * h3 = (TGraph*)obj;
      h3->SetLineColor(color[i]);
      h3->SetMarkerColor(color[i]);
      h3->SetMarkerStyle(style[i]);
      if(draw_) h3->Draw(opt2);
      break;
    }
    default: {
      TString opt1(color[i]==1 ? "" : "same");
      TString opt2 = ("E1X0") + opt1;
      TH1 * h4 = (TH1*)obj;
      h4->SetLineColor(color[i]);
      h4->SetMarkerColor(color[i]);
      h4->SetMarkerStyle(style[i]);
      if(draw_) h4->Draw(opt2);
      break;
    }
    }
  }
}
