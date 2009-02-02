#include <stdio.h>

TList *
makeFileCollection(TString &title="")
{
  ifstream file;
  file.open(title);

  Int_t newline=0;
  char c;
  char line[1000];

  TList * fileCollection = new TList();

  while( file.getline(line,1000)) {
    newline++;
    TFile *f;
    TString name(line);
    f = new TFile(name,"READ");
    fileCollection->Add(f);
  }
  return fileCollection;
}

TList *
makeDirectoryCollection(TList * fileList_ , TString * inDirs, Int_t size)
{
  TList * dirList = new TList();

  TFile * file;
  TDirectory * dir;

  TIter iter(fileList_);
  while ( (file = (TFile *)iter()) ) {
    for (int i = 0; i<size; ++i) {
      TString theDir(inDirs[i]);
      dir = getDirectory(file,theDir);
      dirList->Add(dir);
    }
  }
  return dirList;
}

TList *
makeObjectCollection(TList * dirList_ , TString & histoName_)
{
  TList * objList = new TList();
  TIter iter(dirList_);
  TDirectory * dir;
  while ( (dir = (TDirectory *)iter()) ) {
    TObject * obj = (dir) ? dir->Get(histoName_) : 0;
    objList->Add(obj);
  }
  return objList;
}


TList *
makeFitCollection(TList * objectList_,int fitType = 2) 
{
  TList * fitList = new TList();

  TGraph * fitGraph;
  //pair<double,double> fitPair;
  TObjString * fitString;

  TIter iter(objectList_);
  TObject * obj;
  int i = 0;
  while ( (obj=(TObject*)iter()) )  {
    i++;
    TString iString(i);
    if( obj->IsA()->InheritsFrom("TH2") ) {
      TH2 * h2 = (TH2*)obj;
      fitGraph = fit(h2,fitType);
      fitGraph->SetName(fitGraph->GetName()+iString);
      fitList->Add(fitGraph);
    }
    else if( obj->IsA()->InheritsFrom("TH1") ) {
      TH1 * h1 = (TH1*)obj;
      //fitPair = 
      fit(h1,fitType);
      //fitGraph->SetName(fitGraph->GetName()+iString);
      //fitList->Add(fitGraph);
    }
    else if( obj->IsA()->InheritsFrom("TGraph") ) {
      TGraph * g1 = (TGraph*)obj;
      std::pair<double,double> fitTest = fit(g1,fitType);
      char buffer [50];
      sprintf(buffer,"$%7.4f \\pm %4.2e $",100*fitTest.first,100*fitTest.second);
      //printf("$%7.4f \\pm %4.2e $",100*fitTest.first,100*fitTest.second);

      fitString = new TObjString(buffer);
      fitList->Add(fitString);
    }
  }
  fitList->SetOwner();
  return fitList;
}

TList *
makeTGfromTHCollection(TList * objectlist_)
{
  TList * graphList = new TList();

  TObject * obj;
  TIter iter(objectlist_);
  while ( (obj=(TObject*)iter()) ) {
    if( ! (obj->IsA()->InheritsFrom("TH1") && ! obj->IsA()->InheritsFrom("TH2") ) ) {
      cout << "Not a TH1!" << endl;
      continue;
    }
    TH1* h1 = (TH1*)obj;
    graphList->Add(convertTH1toTGraph(h1));
  }
  graphList->SetOwner();
  return graphList;
}
