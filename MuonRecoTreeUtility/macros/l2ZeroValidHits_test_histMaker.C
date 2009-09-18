#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
//#include "Math/LorentzVector.h"
#include <vector>
#include "CMS1.h"
#include <TStopwatch.h>

//int ScanTree ( TTree* tree) {
int ScanTree ( TTree* tree, char *fileName) {

  // This reads in the tree.  As you might imagine.
  Init(tree);
  TFile *histFile = new TFile(fileName,"RECREATE");
  TDirectory *histDir = histFile->mkdir("eff_hist");
  TStopwatch timer;
  timer.Start();

  // Declare the usual constants
  double pi = 3.141592653;

  histDir->cd();

  int nEntries = tree->GetEntries();
    
  TH1F *troubledL2Parent = new TH1F("troubledL2Parent","MotherBin L2 with no valid hits",13,-0.5,12.5);
  TH1F *troubledL2ParentPt = new TH1F("troubledL2ParentPt","p_{T,sim} for parent, no valid L2 hits",50,0,50);
  TH2F *troubledL2EtaScatter = new TH2F("troubledL2EtaScatter","sim vs L2 #eta, no valid L2 hits",100,-2.4,2.4,100,-2.4,2.4);
  TH1F *troubledL2DR = new TH1F("troubledL2DR","#Delta R(sim, L2), no valid hits",100,0,2);

  //Event Loop
  for( int iEntry = 0; iEntry < nEntries; iEntry++) {
    if (iEntry%1000 == 0) cout << "Event " << iEntry << " time " << timer.RealTime() << " cpu " << timer.CpuTime() << endl;
    timer.Continue();
    tree->GetEntry(iEntry);
    
    try {
      for (int iL2 = 0; iL2 < nL2; iL2++) {
	//		cout << "testing a fucking iterator" << endl;
	for (map<int, int>::const_iterator hits = (*l2NMuHits).begin(); hits != (*l2NMuHits).end(); hits++) {
	  //	  	  cout << "testing if we have no valid hits" << endl;
	  if (hits->first == iL2 && hits->second == 0) { // No valid hits in mu system
	    //	    	    cout << "trying to get the fucking sim index" << endl;
	    int simIndex = RTS_Association_Index((*l2Eta).at(iL2), (*l2Phi).at(iL2), simMuonEta, simMuonPhi, 9.0);
	    double DR = RTS_Association_DeltaR((*l2Eta).at(iL2), (*l2Phi).at(iL2), simMuonEta, simMuonPhi, 9.0);
	    //	    if ((*simMuonMotherBinNumber).at(simIndex) < 0) (*simMuonMotherBinNumber).at(simIndex) = 0;
	    if (simIndex >= 0) {
	      troubledL2Parent->Fill((*simMuonMotherBinNumber).at(simIndex));
	      troubledL2DR->Fill(DR);
	      troubledL2ParentPt->Fill((*simMuonPt).at(simIndex));
	      troubledL2EtaScatter->Fill((*l2Eta).at(iL2),(*simMuonEta).at(simIndex));
	    }
	  }
	}
      }
    }
    catch (...) {
      cout << "Fuck it.  Pressing on" << endl;
    }
  }
  
  histDir->Write("",TObject::kOverwrite);
  
  return 0;
}

void divide_histos_and_errors(TH1F* HIS1, TH1F* HIS2, TH1F* HIS3) {
  // HIS1 = numerator
  // HIS2 = denominator
  // HIS3 = ratio
  Int_t nbins1 = HIS1->GetNbinsX(); // find out nr of bins to loop over.
  for(Int_t ibins1=1 ; ibins1 <= nbins1 ; ibins1 ++) {
    Float_t content_1 = HIS1->GetBinContent(ibins1);
    Float_t content_2 = HIS2->GetBinContent(ibins1);
    Float_t error_1 = HIS1->GetBinError(ibins1);
    Float_t error_2 = HIS2->GetBinError(ibins1);
    
    if (content_2 != 0) {
      HIS3->SetBinContent(ibins1, content_1/content_2);
      HIS3->SetBinError(ibins1, error_1/content_2);
    }
  }
}

void divide_TH2_histos_and_errors(TH2F* HIS1, TH2F* HIS2, TH2F* HIS3) {
  // HIS1 = numerator
  // HIS2 = denominator
  // HIS3 = ratio
  Int_t nbinsX = HIS1->GetNbinsX(); // find out nr of bins to loop over.
  Int_t nbinsY = HIS1->GetNbinsY();
  for(Int_t ibinsX=1; ibinsX <= nbinsX; ibinsX ++) {
    for (Int_t ibinsY=1; ibinsY <= nbinsY; ibinsY ++) {
      Float_t content_1 = HIS1->GetBinContent(ibinsX, ibinsY);
      Float_t content_2 = HIS2->GetBinContent(ibinsX, ibinsY);
      Float_t error_1 = HIS1->GetBinError(ibinsX, ibinsY);
      Float_t error_2 = HIS2->GetBinError(ibinsX, ibinsY);
      
      if (content_2 != 0) {
	HIS3->SetBinContent(ibinsX, ibinsY, content_1/content_2);
	HIS3->SetBinError(ibinsX, ibinsY, error_1/content_2);
      }
      else {
	HIS3->SetBinContent(ibinsX, ibinsY, 0);
	HIS3->SetBinError(ibinsX, ibinsY, 0);
      }
    }      
  }
}

int STR_Association_Index (double simEta, double simPhi, vector<double> *recoEta, vector<double> *recoPhi, double DR_CUT) {
  int index = -1;
  double DeltaR = 999.9;

  //  cout << "testing this method...looping over sim" << endl;
  //  cout << "simEta size" << (*simEta).size() << endl;

  for (int i = 0; i < (*recoEta).size(); i++) {
    //    cout << "eta phi at index" << (*simEta).at(i) <<" "<< (*simPhi).at(i) << endl;
    double temp_deleta = simEta - (*recoEta).at(i);
    double temp_delphi = simPhi - (*recoPhi).at(i);
    double temp_DeltaR = sqrt((temp_deleta * temp_deleta) + (temp_delphi * temp_delphi));
    if (temp_DeltaR < DeltaR  && temp_DeltaR < DR_CUT) {
      index = i;
      DeltaR = temp_DeltaR;
    }
  }
  //  cout << "DeltaR = " << DeltaR << endl;
  //  cout << "returning" << endl;
  return index;
}

int RTS_Association_Index (double recoEta, double recoPhi, vector<double> *simEta, vector<double> *simPhi, double DR_CUT) {
  int index = -1;
  double DeltaR = 999.9;

  //  cout << "testing this method...looping over sim" << endl;
  //  cout << "simEta size" << (*simEta).size() << endl;

  for (int i = 0; i < (*simEta).size(); i++) {
    //    cout << "eta phi at index" << (*simEta).at(i) <<" "<< (*simPhi).at(i) << endl;
    //    cout << "to match" << recoEta << " " << recoPhi << endl;
    double temp_deleta = recoEta - (*simEta).at(i);
    double temp_delphi = recoPhi - (*simPhi).at(i);
    double temp_DeltaR = sqrt((temp_deleta * temp_deleta) + (temp_delphi * temp_delphi));
    if (temp_DeltaR < DeltaR  && temp_DeltaR < DR_CUT) {
      index = i;
      DeltaR = temp_DeltaR;
    }
  }
  //  cout << "DeltaR = " << DeltaR << endl;
  //  cout << "returning " << index << endl;
  return index;
}

double RTS_Association_DeltaR (double recoEta, double recoPhi, vector<double> *simEta, vector<double> *simPhi, double DR_CUT) {
  int index = -1;
  double DeltaR = 999.9;

  //  cout << "testing this method...looping over sim" << endl;
  //  cout << "simEta size" << (*simEta).size() << endl;

  for (int i = 0; i < (*simEta).size(); i++) {
    //    cout << "eta phi at index" << (*simEta).at(i) <<" "<< (*simPhi).at(i) << endl;
    //    cout << "to match" << recoEta << " " << recoPhi << endl;
    double temp_deleta = recoEta - (*simEta).at(i);
    double temp_delphi = recoPhi - (*simPhi).at(i);
    double temp_DeltaR = sqrt((temp_deleta * temp_deleta) + (temp_delphi * temp_delphi));
    if (temp_DeltaR < DeltaR  && temp_DeltaR < DR_CUT) {
      index = i;
      DeltaR = temp_DeltaR;
    }
  }
  //  cout << "DeltaR = " << DeltaR << endl;
  //  cout << "returning " << index << endl;
  return DeltaR;
}

