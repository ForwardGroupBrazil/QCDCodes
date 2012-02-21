#include <iomanip>
void GetTreeSize(TString FileName, TString TreeName)
{
  TFile *inf          = TFile::Open(FileName);
  TTree *tr           = (TTree*)inf->Get(TreeName);
  TObjArray *branches = (TObjArray*)tr->GetListOfBranches();
  int size(0);
  cout.setf(ios::right);
  int N(branches->GetEntries());
  TH1F *hSize = new TH1F("size","size",N,0,N);
  for(int ib=0;ib<N;ib++) {
    TString name(branches->At(ib)->GetName());
    TBranch *br = (TBranch*)tr->GetBranch(name);
    hSize->Fill(name,br->GetZipBytes()/1e+3); 
    size += br->GetZipBytes();
  } 
  cout<<"Total size: "<<size<<endl;
  for(int ib=0;ib<N;ib++) {
    TString name(branches->At(ib)->GetName());
    TBranch *br = (TBranch*)tr->GetBranch(name);
    float percent = TMath::Ceil(1000*float(br->GetZipBytes())/float(size))/10;
    cout<<ib<<setw(20)<<name<<setw(15)<<br->GetZipBytes()<<" "<<percent<<"%"<<endl;
  }

  TCanvas *can = new TCanvas("TreeSize","TreeSize",1000,400);
  hSize->GetXaxis()->SetTitle("Branch Name");
  hSize->GetXaxis()->SetLabelSize(0.04);
  hSize->GetYaxis()->SetTitle("Size (KB)");
  hSize->SetFillColor(kGray);
  hSize->Draw();
}
