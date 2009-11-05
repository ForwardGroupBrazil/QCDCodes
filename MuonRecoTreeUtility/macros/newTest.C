#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
//#include "Math/LorentzVector.h"
#include <vector>
#include "test.h"

#include <sstream>

template <class T>
inline std::string to_string (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}

int ScanTree ( TTree* tree, char *fileName="histo.root") {

  TFile *histFileAlgo = new TFile(fileName,"RECREATE");
  TDirectory *histDirAlgo = histFileAlgo->mkdir("test_hist");

  Init(tree);

  TH1F* h_pt[6];
  TDirectory* dirs[6];
  TH1F* tmp;
  histDirAlgo->cd();
  for (unsigned int ww=0;ww<6;ww++){
    histDirAlgo->cd();
    std::string ss = to_string(ww);

    TString name("dir");
    name+=ss;
    TString name2("l3TestPt");
    //name2+=ss;
    TDirectory *tmpDir = histDirAlgo->mkdir(name.Data());
    dirs[ww] = tmpDir;
    tmpDir->cd();
    tmp = new TH1F(name2,"histogram of L3 p_{T}",100,0,100);
    tmp->SetDirectory(tmpDir);
    h_pt[ww] = tmp;
  }
  histDirAlgo->cd();

  TH1F *l3TestPtAll = new TH1F("l3TestPtAll","histogram of L3 p_{T}",100,0,100);
  TH1F *l3deltaPt = new TH1F("l3deltaPt","histogram of L3 #Delta p_{T}",100,-10.,10.);


  int nEntries = tree->GetEntries();

  //Event Loop
  for( int i = 0; i < nEntries; i++) {
    tree->GetEntry(i);
    //aaa if(nL3 > 0) {
    for (int iL3 = 0; iL3 < nL3; iL3++) {
      /*
	cout << iL3 << " of " << nL3 << " MyBit " << (*l3AssociationMyBit).at(iL3) << endl;
	if((*l3AssociationMyBit).at(iL3) & 1<<0) cout << "0" << endl;
	if((*l3AssociationMyBit).at(iL3) & 1<<1) cout << "1" << endl;
	if((*l3AssociationMyBit).at(iL3) & 1<<2) cout << "2" << endl;
	if((*l3AssociationMyBit).at(iL3) & 1<<3) cout << "3" << endl;
	if((*l3AssociationMyBit).at(iL3) & 1<<4) cout << "4" << endl;
	if((*l3AssociationMyBit).at(iL3) & 1<<5) cout << "5" << endl;
      */
      l3TestPtAll->Fill((*l3Pt).at(iL3));

      for(unsigned int w=0 ; w<6 ; w++){
	int ii = w;
      	if((*l3AssociationMyBit).at(iL3) & 1<<ii ) {
	  h_pt[ii]->Fill((*l3Pt).at(iL3));
	}
      }
      l3deltaPt->Fill((*l2Pt).at(iL3)-(*tkTrackPt).at(iL3));
    }
    //aaa } //cout << "line 83"<< endl;
  }

  for(int j =0; j<6;j++){
    //h_pt[j]->Write("",TObject::kOverwrite);
    dirs[j]->Write("",TObject::kOverwrite);
  }

  histDirAlgo->Write("",TObject::kOverwrite);
  histFileAlgo->Close();

  return 0;
}
