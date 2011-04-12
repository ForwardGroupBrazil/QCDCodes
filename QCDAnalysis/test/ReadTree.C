void ReadTree()
{
  gROOT->ProcessLine("#include <vector>");
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  TFile *inf  = TFile::Open("ProcessedTree_data.root");
  TFile *outf = new TFile("Histo.root","RECREATE");
  TTree *tr = (TTree*)inf->Get("ak5/ProcessedTree");
  //----------------------------------------------------
  TH1F *hPFCorMjj = new TH1F("PFCorMjj","PFCorMjj",500,0,5000);
  //-----------------------------------------------------
  QCDEvent *Event = new QCDEvent();
  TBranch *branch = tr->GetBranch("event");
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
    cout<<"Event #"<<i<<":"<<endl;
    cout<<"PFMass = "<<Event->pfmjjcor(0)<<", JEC up:  "<<Event->pfmjjcor(1)<<", JEC down: "<<Event->pfmjjcor(-1)<<endl;
    cout<<"PFJets: "<<Event->nPFJets()<<endl;
    for(unsigned j=0;j<Event->nPFJets();j++) {
      cout<<j<<" pt = "<<(Event->pfjet(j)).ptCor()<<", cor = "<<(Event->pfjet(j)).cor()<<", unc = "<<(Event->pfjet(j)).unc()<<endl;
    }
    cout<<"CaloJets: "<<Event->nCaloJets()<<endl;
    for(unsigned j=0;j<Event->nCaloJets();j++) {
      cout<<j<<" pt = "<<(Event->calojet(j)).ptCor()<<", cor = "<<(Event->calojet(j)).cor()<<", unc = "<<(Event->calojet(j)).unc()<<endl;
    }
    cout<<"HLT Objects: "<<Event->nHLTObj()<<endl;
    //for(int j=0;j<Event->nHLTObj();j++) {
      //cout<<j<<" pt = "<<(Event->hltobj(j)).pt()<<endl;
    //}
    cout<<"L1 Objects: "<<Event->nL1Obj()<<endl;
    //for(int j=0;j<Event->nL1Obj();j++) {
      //cout<<j<<" pt = "<<(Event->l1obj(j)).pt()<<endl;
    //}
    hPFCorMjj->Fill(Event->pfmjjcor(0));
  }
  //hPFCorMjj->Draw();
  outf->Write();
}
