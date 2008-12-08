#include <stdio.h>

TSortedList *
makeFileCollection(TString &title="")
{
  ifstream file;
  file.open(title);

  Int_t newline=0;
  char c;
  char line[1000];

  TSortedList * fileCollection = new TSortedList();

  while( file.getline(line,1000)) {
    newline++;
    TFile *f;
    TString name(line);
    f = new TFile(name,"READ");
    fileCollection->Add(f);
  }
  return fileCollection;
}

TSortedList *
makeDirectoryCollection(TSortedList * fileList_ , TString * inDirs, Int_t size)
{
  TSortedList * dirList = new TSortedList();

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

TSortedList *
makeObjectCollection(TSortedList * dirList_ , TString & histoName_)
{
  TSortedList * objList = new TSortedList();
  TIter iter(dirList_);
  TDirectory * dir;
  while ( (dir = (TDirectory *)iter()) ) {
    TObject * obj = (dir) ? dir->Get(histoName_) : 0;
    objList->Add(obj);
  }
  return objList;
}


TSortedList *
makeFitCollection(TSortedList * objectList_,int fitType = 2) 
{
  TSortedList * fitList = new TSortedList();

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
