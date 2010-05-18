#include "TH1F.h"

#include <vector>
#include "RecoTree.h"

int ScanTreeBasicsNote ( TTree* tree, char *fileName, bool isData=false, double weight = 1.0) {

  // This reads in the tree.  As you might imagine.
  Init(tree);
  TFile *histFile = new TFile(fileName,"UPDATE");
  TDirectory *histDir = histFile->mkdir("basic_hist");

  histDir->cd();

  TH1F *hist_L3_pt = new TH1F("hist_L3_pt"," p_{T}",20,0.,20.);
  hist_L3_pt->Sumw2();
  TH1F *hist_tk_pt = new TH1F("hist_tk_pt"," p_{T}",20,0.,20.);
  hist_tk_pt->Sumw2();

  TH1F *hist_L3_pt_res = new TH1F("hist_L3_pt_res"," #sigma(#delta p_{T}/p_{T}) [%]",300,-3.,3.);
  TH1F *hist_tk_pt_res = new TH1F("hist_tk_pt_res"," #sigma(#delta p_{T}/p_{T}) [%]",300,-3.,3.);

  TH1F *hist_L3_p = new TH1F("hist_L3_p"," p",20,0.,20.);
  hist_L3_p->Sumw2();
  TH1F *hist_tk_p = new TH1F("hist_tk_p"," p",20,0.,20.);
  hist_tk_p->Sumw2();

  TH1F *hist_L3_eta = new TH1F("hist_L3_eta"," #eta",60,-3.,3.);
  hist_L3_eta->Sumw2();
  TH1F *hist_tk_eta = new TH1F("hist_tk_eta"," #eta",60,-3.,3.);
  hist_tk_eta->Sumw2();

  TH1F *hist_L3_d0 = new TH1F("hist_L3_d0"," p_{T}",20,-1.,1.);
  hist_L3_d0->Sumw2();
  TH1F *hist_tk_d0 = new TH1F("hist_tk_d0"," p_{T}",20,-1.,1.);
  hist_tk_d0->Sumw2();

  TH1F *hist_tk_nHits = new TH1F("hist_tk_nHits","Number of hits",76,-0.5,75.5);
  TH1F *hist_tk_chi2 = new TH1F("hist_tk_chi2","#chi^{2}",100,0.,50.);
  TH1F *hist_tk_nchi2 = new TH1F("hist_tk_nchi2","#chi^{2}_{norm}",100,0.,10.);

  int nEntries = tree->GetEntries();

  //Event Loop
  for( int iEntry = 0; iEntry < nEntries; iEntry++) {
    if (iEntry%1000 == 0) cout << "Event " << iEntry << endl;
    tree->GetEntry(iEntry);
    if(nMu==0) continue;
    
    for(int iMu = 0; iMu < nMu; iMu++) {
      
      bool passL3 = false;
      bool passL2 = false;
      bool passTK = false;
      
      if(    (*muAllGlobalMuons).at(iMu)
	     ) {passL3 = true; passTK = true;}
      
      int myIsAssociated = (*l3IsAssociated).at(iMu);
      int myL3AssociationPdgId = (*l3AssociationPdgId).at(iMu);
      int myL3ParentID = (*l3ParentID).at(iMu);
      int myL3MotherBinNumber = (*l3MotherBinNumber).at(iMu);
      int myL3TPIdBin = GetBinNum((*l3AssociationPdgId).at(iMu));
      int myL3ParentIdBin = GetBinNum((*l3ParentID).at(iMu));
      


      if (passL3) {
	int l3_index = iMu;
	double pt_L3 = (*l3Pt).at(iMu);
	double pt_L3sim = (*l3AssociatedSimMuonPt).at(iMu);
	double p_L3 = (*l3P).at(iMu);
	double d0_L3 = (*l3D0).at(iMu);
	double eta_L3 = (*l3Eta).at(iMu);
	hist_L3_pt ->Fill(pt_L3,!isData?1.0*weight:weight);
	hist_L3_pt_res ->Fill((pt_L3-pt_L3sim)/pt_L3sim,
			      !isData?1.0*weight:weight);
	hist_L3_p ->Fill(p_L3,!isData?1.0*weight:weight);
	hist_L3_d0 ->Fill(d0_L3,!isData?1.0*weight:weight);
	hist_L3_eta->Fill(eta_L3,!isData?1.0*weight:weight);
      }
      
      if (passTK) {
	hist_tk_pt ->Fill((*tkTrackPt).at(iMu),!isData?1.0*weight:weight);
	hist_tk_pt_res ->Fill(((*tkTrackPt).at(iMu)-(*tkTrackAssociatedSimMuonPt).at(iMu))/(*tkTrackAssociatedSimMuonPt).at(iMu),!isData?1.0*weight:weight);
	hist_tk_p ->Fill((*tkTrackP).at(iMu),!isData?1.0*weight:weight);
	hist_tk_eta->Fill((*tkTrackEta).at(iMu),!isData?1.0*weight:weight);
	hist_tk_d0 ->Fill((*tkTrackD0).at(iMu),!isData?1.0*weight:weight);
	hist_tk_nHits->Fill((*tkTrackNHits).at(iMu),!isData?1.0*weight:weight);
	hist_tk_chi2->Fill((*tkTrackChi2).at(iMu),!isData?1.0*weight:weight);
	hist_tk_nchi2->Fill((*tkTrackChi2).at(iMu)/(*tkTrackNdof).at(iMu),!isData?1.0*weight:weight);
      }
      
    }
    
  }
  
  histDir->cd();
  histDir->Write("",TObject::kOverwrite);
  histFile->Close();
  
  return 0;
  
}
