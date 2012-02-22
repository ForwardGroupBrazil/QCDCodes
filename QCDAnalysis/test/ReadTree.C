void ReadTree()
{
  gROOT->ProcessLine("#include <vector>");
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  TFile *inf  = TFile::Open("ProcessedTree_data.root");
  TTree *tr = (TTree*)inf->Get("ak7/ProcessedTree");
  QCDEvent *Event = new QCDEvent();
  TBranch *branch = tr->GetBranch("events");
  branch->SetAddress(&Event);
  
  unsigned NEntries = tr->GetEntries();
  cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
  int    decade = 0;
  
  for(unsigned i=0;i<NEntries;i++) {
    double progress = 10.0*i/(1.0*NEntries);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;          
    tr->GetEntry(i);
    for(unsigned itrig=0;itrig<Event->nTriggers();itrig++) {
      cout<<"Trigger #"<<itrig<<", decision = "<<Event->fired(itrig)<<endl;
      if (Event->fired(itrig) > 0) {
        cout<<"Found "<<Event->nL1Obj(itrig)<<" L1 objects: "<<endl;
        for(unsigned iobj=0;iobj<Event->nL1Obj(itrig);iobj++)
          cout<<"pt = "<<Event->l1obj(itrig,iobj).pt()<<endl;
        cout<<"Found "<<Event->nHLTObj(itrig)<<" HLT objects: "<<endl;
        for(unsigned iobj=0;iobj<Event->nHLTObj(itrig);iobj++)
          cout<<"pt = "<<Event->hltobj(itrig,iobj).pt()<<endl;
      }
    }
  }
}
