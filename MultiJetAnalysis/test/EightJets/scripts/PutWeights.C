void PutWeights(TString FileName, float Weight)
{
  TFile *inf   = TFile::Open(FileName+".root");
  TFile *outf  = TFile::Open(FileName+"_weights.root","RECREATE");
  TTree *trIN  = (TTree*)inf->Get("multijets/events");
  TTree *trOUT = (TTree*)trIN->CloneTree();
  float wt;
  TBranch *brWtMC = trOUT->Branch("wt",&wt,"wt/F");
  for(int i=0;i<trOUT->GetEntries();i++) {
    trOUT->GetEntry(i);
    wt = Weight;
    //mcWt = 1./104.;
    //mcWt = 1./2457.;
    //mcWt = 1./49713.;
    //cout<<weight<<" "<<mcWt<<endl;
    brWtMC->Fill();
  }
  TDirectoryFile *dir = (TDirectoryFile*)outf->mkdir("multijets");
  dir->cd();
  trOUT->Write("events");
  //inf->Close();
}
