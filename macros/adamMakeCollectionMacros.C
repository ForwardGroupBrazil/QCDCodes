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

  TIter iter(objectList_);
  TObject * obj;
  while ( (obj=(TObject*)iter()) )  {
    if( obj->IsA()->InheritsFrom("TH2") ) {
      TH2 * h2 = (TH2*)obj;
      fitGraph = fit(h2,fitType);
      fitList->Add(fitGraph);
    }
  }
  return fitList;
}
