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

  int nPtBins = 19;
  Double_t pt_Edges[20] = {0, 5, 7, 9, 11, 13, 15, 18, 21, 24, 27, 30, 34, 38, 42, 50, 60, 70, 80, 100};
  
  TH1F *l3OverXPtEffNum = new TH1F("l3OverXPtEffNum","l3 associated sim pt",nPtBins,pt_Edges);
  TH1F *l3OverXEtaEffNum = new TH1F("l3OverXEtaEffNum","l3 associated sim #eta",24,-2.4,2.4);
  TH2F *l3OverXPtEtaEffNum = new TH2F("l3OverXPtEtaEffNum","blah",nPtBins,pt_Edges,24,-2.4,2.4);

  TH1F *l3TkOverSimPtEffNum = new TH1F("l3TkOverSimPtEffNum","l3 associated sim pt",nPtBins,pt_Edges);
  TH1F *l3TkOverSimEtaEffNum = new TH1F("l3TkOverSimEtaEffNum","l3 associated sim #eta",24,-2.4,2.4);
  TH2F *l3TkOverSimPtEtaEffNum = new TH2F("l3TkOverSimPtEtaEffNum","blah",nPtBins,pt_Edges,24,-2.4,2.4);

  TH1F *l3OverSimPtEffDenom = new TH1F("l3OverSimPtEffDenom","l3 associated sim pt",nPtBins,pt_Edges);
  TH1F *l3OverSimEtaEffDenom = new TH1F("l3OverSimEtaEffDenom","l3 associated sim #eta",24,-2.4,2.4);
  TH2F *l3OverSimPtEtaEffDenom = new TH2F("l3OverSimPtEtaEffDenom","blah",nPtBins,pt_Edges,24,-2.4,2.4);
  TH1F *l3OverL1PtEffDenom = new TH1F("l3OverL1PtEffDenom","l3 associated sim pt",nPtBins,pt_Edges);
  TH1F *l3OverL1EtaEffDenom = new TH1F("l3OverL1EtaEffDenom","l3 associated sim #eta",24,-2.4,2.4);
  TH2F *l3OverL1PtEtaEffDenom = new TH2F("l3OverL1PtEtaEffDenom","blah",nPtBins,pt_Edges,24,-2.4,2.4);
  TH1F *l3OverL2PtEffDenom = new TH1F("l3OverL2PtEffDenom","l3 associated sim pt",nPtBins,pt_Edges);
  TH1F *l3OverL2EtaEffDenom = new TH1F("l3OverL2EtaEffDenom","l3 associated sim #eta",24,-2.4,2.4);
  TH2F *l3OverL2PtEtaEffDenom = new TH2F("l3OverL2PtEtaEffDenom","blah",nPtBins,pt_Edges,24,-2.4,2.4);

  TH1F *l3OverSimPtEff = new TH1F("l3OverSimPtEff","l3/sim eff vs. assoc. sim pt",nPtBins,pt_Edges);
  TH1F *l3OverSimEtaEff = new TH1F("l3OverSimEtaEff","l3/sim eff vs. assoc. sim #eta",24,-2.4,2.4);
  TH2F *l3OverSimPtEtaEff = new TH2F("l3OverSimPtEtaEff","l3/sim eff map: p_{T} vs #eta",nPtBins,pt_Edges,24,-2.4,2.4);
  TH1F *l3TkOverSimPtEff = new TH1F("l3TkOverSimPtEff","l3/sim eff vs. assoc. sim pt",nPtBins,pt_Edges);
  TH1F *l3TkOverSimEtaEff = new TH1F("l3TkOverSimEtaEff","l3/sim eff vs. assoc. sim #eta",24,-2.4,2.4);
  TH2F *l3TkOverSimPtEtaEff = new TH2F("l3TkOverSimPtEtaEff","l3/sim eff map: p_{T} vs #eta",nPtBins,pt_Edges,24,-2.4,2.4); 
  TH1F *l3OverL1PtEff = new TH1F("l3OverL1PtEff","l3/l1 eff vs. assoc. sim pt",nPtBins,pt_Edges);
  TH1F *l3OverL1EtaEff = new TH1F("l3OverL1EtaEff","l3/l1 eff vs. assoc. sim #eta",24,-2.4,2.4);
  TH2F *l3OverL1PtEtaEff = new TH2F("l3OverL1PtEtaEff","l3/l1 eff map: p_{T} vs #eta",nPtBins,pt_Edges,24,-2.4,2.4);
  TH1F *l3OverL2PtEff = new TH1F("l3OverL2PtEff","l3/l2 eff vs. assoc. sim pt",nPtBins,pt_Edges);
  TH1F *l3OverL2EtaEff = new TH1F("l3OverL2EtaEff","l3/l2 eff vs. assoc. sim #eta",24,-2.4,2.4);
  TH2F *l3OverL2PtEtaEff = new TH2F("l3OverL2PtEtaEff","l3/l2 eff map: p_{T} vs #eta",nPtBins,pt_Edges,24,-2.4,2.4);
  
  

  TH1F *muHLTTime = new TH1F("muHLTTime","time taken for Muon HLT",150,0,0.3);
  TH1F *muHLTTime_L2Filtered = new TH1F("muHLTTime_L2Filtered","time taken for Muon HLT, L2 filtering applied",150,0,0.3);

  TH1F *l3SeedTime = new TH1F("time_l3Seed","time taken for hltL3TrajectorySeed",100,0,0.05);
  TH1F *l3SeedTimeIO = new TH1F("time_l3SeedIO","time taken for hltL3TrajectorySeedIO",100,0,0.05);
  TH1F *l3TkCandTime = new TH1F("time_l3TkCand","time taken for hltL3TrackCandidateFromL2",100,0,0.1);
  TH1F *l3TkCandTimeIO = new TH1F("time_l3TkCandIO","time taken for hltL3TrackCandidateFromL2IO",100,0,0.1);
  TH1F *l3TkTime = new TH1F("time_l3Tk","time taken for hltL3TracksFromL2OI",100,0,0.05);
  TH1F *l3TkTimeIO = new TH1F("time_l3TkIO","time taken for hltL3TracksFromL2IO",100,0,0.05);
  TH1F *l3TkTimeMerge = new TH1F("time_l3TkMerge","time taken for hltL3TracksFromL2 (merge)",100,0,0.05);
  TH1F *l3MuTime = new TH1F("time_l3Mu","time taken for hltL3Muons",100,0,0.1);
  TH1F *l3MuCandTime = new TH1F("time_l3MuCand","time taken for hltL3MuonCandidates",100,0,0.05);



  TH1F *l2SeedTime = new TH1F("l2SeedTime","time taken for hltL2MuonSeeds",100,0,0.05);
  TH1F *l2CandTime = new TH1F("l2CandTime","time taken for hltL2MuonCandidates",100,0,0.1);
  TH1F *l2MuTime = new TH1F("l2MuTime","time taken for hltL2Muons",100,0,0.1);  
  TH1F *nCands = new TH1F("nCands","number of Candidates from L2",16,-0.5,15.5);
  TH1F *nSeeds = new TH1F("nSeeds","number of L3 seeds per Event",41,-0.5,40.5);
  TH1F *nSeedsPerL2 = new TH1F("nSeedsPerL2","number of L3 seeds per L2",41,-0.5,40.5);
  TH1F *nTracks = new TH1F("nTracks","number of L3 tracks within search window",16,-0.5,15.5);

  TH1F *l2_qoverp_pull = new TH1F("l2_qoverp_pull","pull of L2 qoverp vs sim",100,-10,10);
  TH1F *l2_lambda_pull = new TH1F("l2_lambda_pull","pull of L2 lambda vs sim",100,-10,10);
  TH1F *l2_phi_pull = new TH1F("l2_phi_pull","pull of L2 phi vs sim",100,-10,10);
  TH1F *l2_dxy_pull = new TH1F("l2_dxy_pull","pull of L2 dxy vs sim",100,-1,1);
  TH1F *l2_dsz_pull = new TH1F("l2_dsz_pull","pull of L2 dsz vs sim",100,-10,10);

  TH2F *l2_qoverp_pull_vsPt = new TH2F("l2_qoverp_pull_vsPt","pull of L2 qoverp vs p_{T}",nPtBins,pt_Edges,100,-10,10);
  TH2F *l2_lambda_pull_vsPt = new TH2F("l2_lambda_pull_vsPt","pull of L2 lambda vs p_{T}",nPtBins,pt_Edges,100,-10,10);
  TH2F *l2_phi_pull_vsPt = new TH2F("l2_phi_pull_vsPt","pull of L2 phi vs p_{T}",nPtBins,pt_Edges,100,-10,10);
  TH2F *l2_dxy_pull_vsPt = new TH2F("l2_dxy_pull_vsPt","pull of L2 dxy vs p_{T}",nPtBins,pt_Edges,100,-1,1);
  TH2F *l2_dsz_pull_vsPt = new TH2F("l2_dsz_pull_vsPt","pull of L2 dsz vs p_{T}",nPtBins,pt_Edges,100,-10,10);

  TH2F *l2_qoverp_pull_vsEta = new TH2F("l2_qoverp_pull_vsEta","pull of L2 qoverp vs #eta",24,-2.4,2.4,100,-10,10);
  TH2F *l2_lambda_pull_vsEta = new TH2F("l2_lambda_pull_vsEta","pull of L2 lambda vs #eta",24,-2.4,2.4,100,-10,10);
  TH2F *l2_phi_pull_vsEta = new TH2F("l2_phi_pull_vsEta","pull of L2 phi vs #eta",24,-2.4,2.4,100,-10,10);
  TH2F *l2_dxy_pull_vsEta = new TH2F("l2_dxy_pull_vsEta","pull of L2 dxy vs #eta",24,-2.4,2.4,100,-1,1);
  TH2F *l2_dsz_pull_vsEta = new TH2F("l2_dsz_pull_vsEta","pull of L2 dsz vs #eta",24,-2.4,2.4,100,-10,10);

  TH1F *h_nL3TracksFromL2 = new TH1F("nL3TracksFromL2","Number of L3 Tk Tracks from L2",21,-0.5,20.5);
  TH1F *h_nL3SeedsFromL2 = new TH1F("nL3SeedsFromL2","Number of L3 Seeds from L2",21,-0.5,20.5);
  TH1F *h_nL3 = new TH1F("nL3","Number of L3",21,-0.5,20.5);
  TH1F *h_L3nHits = new TH1F("L3nHits","Number of L3 Hits",76,-0.5,75.5);
  TH1F *h_L3normChi2 = new TH1F("L3normChi2","Normalized #Chi^{2}",200,0.,200.);

  int nEntries = tree->GetEntries();
    
  //Event Loop
  for( int iEntry = 0; iEntry < nEntries; iEntry++) {
    if (iEntry%10000 == 0) cout << "Event " << iEntry << " time " << timer.RealTime() << " cpu " << timer.CpuTime() << endl;
    timer.Continue();
    tree->GetEntry(iEntry);
    
    double max_simPt = -999;

    // Fill in the pure reco crap
    muHLTTime->Fill(totalMuonHLTTime);
    nCands->Fill(nL3Cands);
    nSeeds->Fill(nL3Seeds);
    nTracks->Fill(nL3TracksFromL2);
    for (int i = 0; i < nL2; i++) {
      nSeedsPerL2->Fill((*l2NSeeds).at(i));
    }

    map<std::string,double>::iterator fuckIt;

    try {
      //      cout << "Trying times" << endl;
      for (fuckIt = (*muonL3RecModuleTimes).begin(); fuckIt != (*muonL3RecModuleTimes).end(); fuckIt++) {
	if(iEntry == 10) cout << "Time Iterator a: " << fuckIt->first << endl;
	if (fuckIt->first == "hltL3TrajectorySeed") l3SeedTime->Fill(fuckIt->second);
	if (fuckIt->first == "hltL3TrajSeedIOHit") l3SeedTimeIO->Fill(fuckIt->second);
	if (fuckIt->first == "hltL3TrackCandidateFromL2") l3TkCandTime->Fill(fuckIt->second);
	if (fuckIt->first == "hltL3TrackCandidateFromL2IO") l3TkCandTimeIO->Fill(fuckIt->second);
	if (fuckIt->first == "hltL3TkTracksFromL2OI") l3TkTime->Fill(fuckIt->second);
	if (fuckIt->first == "hltL3TkTracksFromL2IO") l3TkTimeIO->Fill(fuckIt->second);
	if (fuckIt->first == "hltL3TkTracksFromL2") l3TkTimeMerge->Fill(fuckIt->second);
	if (fuckIt->first == "hltL3Muons") l3MuTime->Fill(fuckIt->second);
	if (fuckIt->first == "hltL3MuonCandidates") l3MuCandTime->Fill(fuckIt->second);
      }
      for (fuckIt = (*muonL2RecModuleTimes).begin(); fuckIt != (*muonL2RecModuleTimes).end(); fuckIt++) {
	if (fuckIt->first == "hltL2Muons") l2MuTime->Fill(fuckIt->second);
	if (fuckIt->first == "hltL2MuonSeeds") l2SeedTime->Fill(fuckIt->second);
	if (fuckIt->first == "hltL2MuonCandidates") l2CandTime->Fill(fuckIt->second);
      }
    }
    catch (...) {
      cout << "Fuck it.  Pressing on" << endl;
    }

    

    try { // Try 2
      //      cout << "trying passL1" << endl;
      bool passL1 = false;
      for (int iL1 = 0; iL1 < nL1; iL1++) {
	if ((*l1Pt).at(iL1) + 0.01 > 0 && (*l1Quality).at(iL1) >= 4) {
	  passL1 = true;
	}
      }

      //      cout << "trying passL2" << endl;
      bool passL2 = false;
      for (int iL2 = 0; iL2 < nL2; iL2++) {
	if (passL1) {
	  if ((*l2Pt).at(iL2) > 3 && (*l2Eta).at(iL2) > -2.5 && (*l2Eta).at(iL2) < 2.5) {
	    passL2 = true;
	  }
	}
      }

      if (passL2) muHLTTime_L2Filtered->Fill(totalMuonHLTTime);

      //Adam start
      h_nL3TracksFromL2->Fill(nL3TracksFromL2);
      h_nL3SeedsFromL2->Fill(nL3Seeds); 
      h_nL3->Fill(nL3); 
      for (int i = 0; i < nL3; i++) {
	h_L3nHits->Fill((*l3NHits).at(i));
	h_L3normChi2->Fill((*l3Chi2).at(i)/(*l3Ndof).at(i));
      }
      //Adam end


      // As a temporary way forward, we're doing l2 and l3 reco->sim while I work out what's wrong with the &(*^&*%$# STR association

      if (passL2) {

	for (int i = 0; i < nL2; i++) {
	  if ((*l2IsAssociated).at(i)) {
	    l3OverL2PtEffDenom->Fill((*l2AssociatedSimMuonPt).at(i));
	    l3OverL2EtaEffDenom->Fill((*l2AssociatedSimMuonEta).at(i));
	    l3OverL2PtEtaEffDenom->Fill((*l2AssociatedSimMuonPt).at(i),(*l2AssociatedSimMuonEta).at(i));
	  }
	}
	
	for (int i = 0; i < nL3; i++) {
	  if ((*l3IsAssociated).at(i)) {
	    int l2SeedIndex = (*indexL2SeedingL3).at(i);
	    if ((*l2IsAssociated).at(l2SeedIndex)) {
	      l3OverXPtEffNum->Fill((*l3AssociatedSimMuonPt).at(i));
	      l3OverXEtaEffNum->Fill((*l3AssociatedSimMuonEta).at(i));
	      l3OverXPtEtaEffNum->Fill((*l3AssociatedSimMuonPt).at(i),(*l3AssociatedSimMuonEta).at(i));
	    }
	  }
	}
	
	
	//      cout << "doing efficiencies, resolutions, and shit." << endl;
	//      cout << "how big is our fucking error matrix? " << (*muonErrorMatrix).size() << endl;
	for (int i = 0; i < nSimMuon; i++) {
	  if ((*simMuonPt).at(i) > 3) {
	    //	cout << "filling l3oversim" << endl;
	    l3OverSimPtEffDenom->Fill((*simMuonPt).at(i));
	    l3OverSimEtaEffDenom->Fill((*simMuonEta).at(i));
	    l3OverSimPtEtaEffDenom->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));
	    
	    //	cout << "trying to associate sim to L1 by DR" << endl;
	    int l1Index = STR_L1_Association_Index((*simMuonEta).at(i),(*simMuonPhi).at(i),l1Eta,l1Phi);
	    if (l1Index >= 0) {
	      // cout << "successful.  filling" << endl;
	      l3OverL1PtEffDenom->Fill((*simMuonPt).at(i));
	      l3OverL1EtaEffDenom->Fill((*simMuonEta).at(i));
	      l3OverL1PtEtaEffDenom->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));
	      //	  cout << "erasing a few thing from L1 to prevent double-counting" << endl;
	      (*l1Pt).erase((*l1Pt).begin() + l1Index);
	      (*l1Eta).erase((*l1Eta).begin() + l1Index);
	      (*l1Phi).erase((*l1Phi).begin() + l1Index);
	    } // l1 simtoreco
	    /*	  
	      if ((*simToL2Associated).at(i)) {
	      int index = (*simToL2RecoIndex).at(i);
	      l3OverL2PtEffDenom->Fill((*simMuonPt).at(i));
	      l3OverL2EtaEffDenom->Fill((*simMuonEta).at(i));
	      l3OverL2PtEtaEffDenom->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));
	      
	      // let's test this pull shit out
	    /*	    map<int,vector<double> >::const_iterator yourMom = (*muonErrorMatrix).begin();
	    for ( ; yourMom != (*muonErrorMatrix).end(); ++yourMom) {
	    if (yourMom->first == index) {
	    double tmp_qoverp_pull = ((*l2Qoverp).at(index) - (*simMuonQoverP).at(i))/sqrt(yourMom->second[0]);
	    double tmp_lambda_pull = ((*l2Lambda).at(index) - (*simMuonLambda).at(i))/sqrt(yourMom->second[6]);
	    double tmp_phi_pull = ((*l2Phi).at(index) - (*simMuonPhi).at(i))/sqrt(yourMom->second[12]);
	    double tmp_dxy_pull = ((*l2Dxy).at(index) - (*simMuonDxy).at(i))/sqrt(yourMom->second[18]);
	    double tmp_dsz_pull = ((*l2Dsz).at(index) - (*simMuonDsz).at(i))/sqrt(yourMom->second[24]);
	    
	    l2_qoverp_pull->Fill(tmp_qoverp_pull);
	    l2_lambda_pull->Fill(tmp_lambda_pull);
	    l2_phi_pull->Fill(tmp_phi_pull);
	    l2_dxy_pull->Fill(tmp_dxy_pull);
	    l2_dsz_pull->Fill(tmp_dsz_pull);
	    
	    l2_qoverp_pull_vsPt->Fill((*simMuonPt).at(i),tmp_qoverp_pull);
	    l2_lambda_pull_vsPt->Fill((*simMuonPt).at(i),tmp_lambda_pull);
	    l2_phi_pull_vsPt->Fill((*simMuonPt).at(i),tmp_phi_pull);
	    l2_dxy_pull_vsPt->Fill((*simMuonPt).at(i),tmp_dxy_pull);
	    l2_dsz_pull_vsPt->Fill((*simMuonPt).at(i),tmp_dsz_pull);
	    
	    l2_qoverp_pull_vsEta->Fill((*simMuonEta).at(i),tmp_qoverp_pull);
	    l2_lambda_pull_vsEta->Fill((*simMuonEta).at(i),tmp_lambda_pull);
	    l2_phi_pull_vsEta->Fill((*simMuonEta).at(i),tmp_phi_pull);
	    l2_dxy_pull_vsEta->Fill((*simMuonEta).at(i),tmp_dxy_pull);
	    l2_dsz_pull_vsEta->Fill((*simMuonEta).at(i),tmp_dsz_pull);
	    }
	    }
	    
	    } // l2 STR association
	      // L3 STR association
	      if ((*simToL3Associated).at(i) && (*simToL2Associated).at(i)) {
	      l3OverXPtEffNum->Fill((*simMuonPt).at(i));
	      l3OverXEtaEffNum->Fill((*simMuonEta).at(i));
	      l3OverXPtEtaEffNum->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));
	      }  
	    */
	    if ((*simToTkAssociated).at(i)) {
	      l3TkOverSimPtEffNum->Fill((*simMuonPt).at(i));
	      l3TkOverSimEtaEffNum->Fill((*simMuonEta).at(i));
	      l3TkOverSimPtEtaEffNum->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));
	    }
	  } //end simMuonPt > 3
	} //end loop over SimMuons
      } //end passL2
    } //end Try 2
    catch (...) {
      cout << "Fuck it.  Pressing on" << endl;
    }
  }
  
  divide_histos_and_errors(l3OverXPtEffNum, l3OverSimPtEffDenom, l3OverSimPtEff);
  divide_histos_and_errors(l3OverXEtaEffNum, l3OverSimEtaEffDenom, l3OverSimEtaEff);  
  divide_TH2_histos_and_errors(l3OverXPtEtaEffNum, l3OverSimPtEtaEffDenom, l3OverSimPtEtaEff);  
  
  divide_histos_and_errors(l3OverXPtEffNum, l3OverL1PtEffDenom, l3OverL1PtEff);
  divide_histos_and_errors(l3OverXEtaEffNum, l3OverL1EtaEffDenom, l3OverL1EtaEff);  
  divide_TH2_histos_and_errors(l3OverXPtEtaEffNum, l3OverL1PtEtaEffDenom, l3OverL1PtEtaEff);  
  
  divide_histos_and_errors(l3OverXPtEffNum, l3OverL2PtEffDenom, l3OverL2PtEff);
  divide_histos_and_errors(l3OverXEtaEffNum, l3OverL2EtaEffDenom, l3OverL2EtaEff);  
  divide_TH2_histos_and_errors(l3OverXPtEtaEffNum, l3OverL2PtEtaEffDenom, l3OverL2PtEtaEff);  

  divide_histos_and_errors(l3TkOverSimPtEffNum, l3OverSimPtEffDenom, l3TkOverSimPtEff);
  divide_histos_and_errors(l3TkOverSimEtaEffNum, l3OverSimEtaEffDenom, l3TkOverSimEtaEff);  
  divide_TH2_histos_and_errors(l3TkOverSimPtEtaEffNum, l3OverSimPtEtaEffDenom, l3TkOverSimPtEtaEff);  

  histDir->Write("",TObject::kOverwrite);
  histFile->Close();  
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

int STR_L1_Association_Index (double simEta, double simPhi, vector<double> *l1Eta, vector<double> *l1Phi) {
  int index = -1;
  double DeltaR = 999.9;

  //  cout << "testing this method...looping over sim" << endl;
  //  cout << "simEta size" << (*simEta).size() << endl;

  for (int i = 0; i < (*l1Eta).size(); i++) {
    //    cout << "eta phi at index" << (*simEta).at(i) <<" "<< (*simPhi).at(i) << endl;
    double temp_deleta = simEta - (*l1Eta).at(i);
    double temp_delphi = simPhi - (*l1Phi).at(i);
    double temp_DeltaR = sqrt((temp_deleta * temp_deleta) + (temp_delphi * temp_delphi));
    if (temp_DeltaR < DeltaR  && temp_DeltaR < 1) {
      index = i;
      DeltaR = temp_DeltaR;
    }
  }
  //  cout << "DeltaR = " << DeltaR << endl;
  //  cout << "returning" << endl;
  return index;
}
