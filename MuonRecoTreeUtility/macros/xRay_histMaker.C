#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
//#include "Math/LorentzVector.h"
#include <vector>
#include "RecoTree.h"

//int ScanTree ( TTree* tree) {
int ScanTreePt ( TTree* tree, char *fileName, bool isData=false,double weight = 1.0) {

  // This reads in the tree.  As you might imagine.
  Init(tree);
  TFile *histFile = new TFile(fileName,"UPDATE");
  TDirectory *histDir = histFile->mkdir("PtRate");

  histDir->cd();

  // Time for some histograms and stacks

  double ptLowerLimit = 0;
  double ptUpperLimit = 20;
  double nBins = 20;

  TH1F *l3PtRate_motherBin_1 = new TH1F("l3PtRate_motherBin_1"," #pi^{+/-}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_2 = new TH1F("l3PtRate_motherBin_2"," K",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_3 = new TH1F("l3PtRate_motherBin_3"," D",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_4 = new TH1F("l3PtRate_motherBin_4"," B",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_5 = new TH1F("l3PtRate_motherBin_5"," #Lambda_{b}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_6 = new TH1F("l3PtRate_motherBin_6"," J/#Psi",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_10 = new TH1F("l3PtRate_motherBin_10"," #tau^{+/-}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_11 = new TH1F("l3PtRate_motherBin_11"," #mu",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_12 = new TH1F("l3PtRate_motherBin_12"," b/c",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_13 = new TH1F("l3PtRate_motherBin_13"," other",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_14 = new TH1F("l3PtRate_motherBin_14"," non-associated",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_15 = new TH1F("l3PtRate_motherBin_15"," Data",nBins,ptLowerLimit,ptUpperLimit);

  TH1F *l3PtRate_motherBin_1_cd = new TH1F("l3PtRate_motherBin_1_cd"," #pi^{+/-}",nBins,ptLowerLimit,ptUpperLimit);  
  TH1F *l3PtRate_motherBin_2_cd = new TH1F("l3PtRate_motherBin_2_cd"," K",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_3_cd = new TH1F("l3PtRate_motherBin_3_cd"," D",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_4_cd = new TH1F("l3PtRate_motherBin_4_cd"," B",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_5_cd = new TH1F("l3PtRate_motherBin_5_cd"," #Lambda_{b}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_6_cd = new TH1F("l3PtRate_motherBin_6_cd"," J/#Psi",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_10_cd = new TH1F("l3PtRate_motherBin_10_cd"," #tau^{+/-}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_12_cd = new TH1F("l3PtRate_motherBin_12_cd"," b/c",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_11_cd = new TH1F("l3PtRate_motherBin_11_cd"," #mu",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_13_cd = new TH1F("l3PtRate_motherBin_13_cd"," other",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_14_cd = new TH1F("l3PtRate_motherBin_14_cd"," non-associated",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l3PtRate_motherBin_15_cd = new TH1F("l3PtRate_motherBin_15_cd"," Data",nBins,ptLowerLimit,ptUpperLimit);



  THStack *l3PtRate = new THStack("l3PtRate","L3 rate as f(p_{T,L3})");
  THStack *l3PtRate_cd = new THStack("l3PtRate_cd","L3 rate as f(p_{T,L3})");
  //  THStack *l3PreIsoPtRate = new THStack("l3PreIsoPtRate","L3PreIso rate as f(p_{T,L3})");
  //  THStack *l3PreIsoPtRate_cd = new THStack("l3PreIsoPtRate_cd","L3PreIso rate as f(p_{T,L3})");
  //  THStack *l3IsoPtRate = new THStack("l3IsoPtRate","L3Iso rate as f(p_{T,L3})");
  //  THStack *l3IsoPtRate_cd = new THStack("l3IsoPtRate_cd","L3Iso rate as f(p_{T,L3})");

  TH1F *l2PtRate_motherBin_1 = new TH1F("l2PtRate_motherBin_1"," #pi^{+/-}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_2 = new TH1F("l2PtRate_motherBin_2"," K",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_3 = new TH1F("l2PtRate_motherBin_3"," D",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_4 = new TH1F("l2PtRate_motherBin_4"," B",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_5 = new TH1F("l2PtRate_motherBin_5"," #Lambda_{b}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_6 = new TH1F("l2PtRate_motherBin_6"," J/#Psi",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_10 = new TH1F("l2PtRate_motherBin_10"," #tau^{+/-}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_11 = new TH1F("l2PtRate_motherBin_11"," #mu",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_12 = new TH1F("l2PtRate_motherBin_12"," b/c",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_13 = new TH1F("l2PtRate_motherBin_13"," other",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_14 = new TH1F("l2PtRate_motherBin_14"," non-associated",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_15 = new TH1F("l2PtRate_motherBin_15"," Data",nBins,ptLowerLimit,ptUpperLimit);

  TH1F *l2PtRate_motherBin_1_cd = new TH1F("l2PtRate_motherBin_1_cd"," #pi^{+/-}",nBins,ptLowerLimit,ptUpperLimit);  
  TH1F *l2PtRate_motherBin_2_cd = new TH1F("l2PtRate_motherBin_2_cd"," K",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_3_cd = new TH1F("l2PtRate_motherBin_3_cd"," D",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_4_cd = new TH1F("l2PtRate_motherBin_4_cd"," B",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_5_cd = new TH1F("l2PtRate_motherBin_5_cd"," #Lambda_{b}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_6_cd = new TH1F("l2PtRate_motherBin_6_cd"," J/#Psi",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_10_cd = new TH1F("l2PtRate_motherBin_10_cd"," #tau^{+/-}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_11_cd = new TH1F("l2PtRate_motherBin_11_cd"," #mu",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_12_cd = new TH1F("l2PtRate_motherBin_12_cd"," b/c",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_13_cd = new TH1F("l2PtRate_motherBin_13_cd"," other",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_14_cd = new TH1F("l2PtRate_motherBin_14_cd"," non-associated",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *l2PtRate_motherBin_15_cd = new TH1F("l2PtRate_motherBin_15_cd"," Data",nBins,ptLowerLimit,ptUpperLimit);



  THStack *l2PtRate = new THStack("l2PtRate","L2 rate as f(p_{T,L2})");
  THStack *l2PtRate_cd = new THStack("l2PtRate_cd","L2 rate as f(p_{T,L2})");
  //  THStack *l2PreIsoPtRate = new THStack("l2PreIsoPtRate","L2PreIso rate as f(p_{T,L2})");
  //  THStack *l2PreIsoPtRate_cd = new THStack("l2PreIsoPtRate_cd","L2PreIso rate as f(p_{T,L2})");
  //  THStack *l2IsoPtRate = new THStack("l2IsoPtRate","L2Iso rate as f(p_{T,L2})");
  //  THStack *l2IsoPtRate_cd = new THStack("l2IsoPtRate_cd","L2Iso rate as f(p_{T,L2})");

  TH1F *tkTrackPtRate_motherBin_1 = new TH1F("tkTrackPtRate_motherBin_1"," #pi^{+/-}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_2 = new TH1F("tkTrackPtRate_motherBin_2"," K",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_3 = new TH1F("tkTrackPtRate_motherBin_3"," D",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_4 = new TH1F("tkTrackPtRate_motherBin_4"," B",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_5 = new TH1F("tkTrackPtRate_motherBin_5"," #Lambda_{b}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_6 = new TH1F("tkTrackPtRate_motherBin_6"," J/#Psi",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_10 = new TH1F("tkTrackPtRate_motherBin_10"," #tau^{+/-}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_11 = new TH1F("tkTrackPtRate_motherBin_11"," #mu",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_12 = new TH1F("tkTrackPtRate_motherBin_12"," b/c",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_13 = new TH1F("tkTrackPtRate_motherBin_13"," other",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_14 = new TH1F("tkTrackPtRate_motherBin_14"," non-associated",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_15 = new TH1F("tkTrackPtRate_motherBin_15"," Data",nBins,ptLowerLimit,ptUpperLimit);

  TH1F *tkTrackPtRate_motherBin_1_cd = new TH1F("tkTrackPtRate_motherBin_1_cd"," #pi^{+/-}",nBins,ptLowerLimit,ptUpperLimit);  
  TH1F *tkTrackPtRate_motherBin_2_cd = new TH1F("tkTrackPtRate_motherBin_2_cd"," K",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_3_cd = new TH1F("tkTrackPtRate_motherBin_3_cd"," D",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_4_cd = new TH1F("tkTrackPtRate_motherBin_4_cd"," B",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_5_cd = new TH1F("tkTrackPtRate_motherBin_5_cd"," #Lambda_{b}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_6_cd = new TH1F("tkTrackPtRate_motherBin_6_cd"," J/#Psi",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_10_cd = new TH1F("tkTrackPtRate_motherBin_10_cd"," #tau^{+/-}",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_11_cd = new TH1F("tkTrackPtRate_motherBin_11_cd"," #mu",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_12_cd = new TH1F("tkTrackPtRate_motherBin_12_cd"," b/c",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_13_cd = new TH1F("tkTrackPtRate_motherBin_13_cd"," other",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_14_cd = new TH1F("tkTrackPtRate_motherBin_14_cd"," non-associated",nBins,ptLowerLimit,ptUpperLimit);
  TH1F *tkTrackPtRate_motherBin_15_cd = new TH1F("tkTrackPtRate_motherBin_15_cd"," Data",nBins,ptLowerLimit,ptUpperLimit);



  THStack *tkTrackPtRate = new THStack("tkTrackPtRate","tkTrack rate as f(p_{T,tkTrack})");
  THStack *tkTrackPtRate_cd = new THStack("tkTrackPtRate_cd","tkTrack rate as f(p_{T,tkTrack})");
  //  THStack *tkTrackPreIsoPtRate = new THStack("tkTrackPreIsoPtRate","tkTrackPreIso rate as f(p_{T,tkTrack})");
  //  THStack *tkTrackPreIsoPtRate_cd = new THStack("tkTrackPreIsoPtRate_cd","tkTrackPreIso rate as f(p_{T,tkTrack})");
  //  THStack *tkTrackIsoPtRate = new THStack("tkTrackIsoPtRate","tkTrackIso rate as f(p_{T,tkTrack})");
  //  THStack *tkTrackIsoPtRate_cd = new THStack("tkTrackIsoPtRate_cd","tkTrackIso rate as f(p_{T,tkTrack})");

  // There are many comments to be made about ROOT.  Most of them Rated-R or Rated-X.
  // Here, I'll limit myself to saying this is necessary for the TStacks to go into the TDirectory.
  gDirectory->Append(l2PtRate);
  gDirectory->Append(l3PtRate);
  gDirectory->Append(tkTrackPtRate);

  gDirectory->Append(l2PtRate_cd);
  gDirectory->Append(l3PtRate_cd);
  gDirectory->Append(tkTrackPtRate_cd);

  //  gDirectory->Append(l2IsoPtRate);
  //  gDirectory->Append(l3PreIsoPtRate);
  //  gDirectory->Append(l3IsoPtRate);

  //  gDirectory->Append(l2IsoPtRate_cd);
  //  gDirectory->Append(l3PreIsoPtRate_cd);
  //  gDirectory->Append(l3IsoPtRate_cd);

  int l2OverFlow[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int l3OverFlow[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int tkTrackOverFlow[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  //int l2IsoOverFlow[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
  //int l3PreIsoOverFlow[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
  //int l3IsoOverFlow[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};

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
      
      //if ((*muAllGlobalMuons).at(iMu)) passL3 = true;

      //if ((*muAllStandAloneMuons).at(iMu)) passL2 = true;

      //if ((*muAllGlobalMuons).at(iMu)) passTK = true;

      if (   (*muAllGlobalMuons).at(iMu)
       ) {passL3 = true; passTK = true;}

      //if(passL3) cout << "line 160 " << (*l3Pt).at(iMu) << endl;

      if (passL3) {
	//double pt_L3 = findMaxPt(l3Pt,l3Eta,l3D0);
	//int l3_index = findIndexOfMaxPt(l3Pt,l3Eta,l3D0);
	int l3_index = iMu;
	//if(!(*muAllGlobalMuons).at(iMu)) break;
	double pt_L3 = (*l3Pt).at(iMu);
	if(!isData){
	  if ((*l3IsAssociated).at(l3_index) == 1 && (*l3AssociationVar).at(l3_index) > -0.1 && fabs((*l3AssociationPdgId).at(l3_index))==13) {
	  switch ( (*l3MotherBinNumber).at(l3_index) ) {
	  case 1 :
	    l3PtRate_motherBin_1->Fill(pt_L3);
	    if (pt_L3 > ptUpperLimit) l3OverFlow[0] ++;
	    break;
	  case 2 :
	    l3PtRate_motherBin_2->Fill(pt_L3);
	    if (pt_L3 > ptUpperLimit) l3OverFlow[1] ++;
	    break;
	  case 3 :
	    l3PtRate_motherBin_3->Fill(pt_L3);
	    if (pt_L3 > ptUpperLimit) l3OverFlow[2] ++;
	    break;
	  case 4 :
	    l3PtRate_motherBin_4->Fill(pt_L3);
	    if (pt_L3 > ptUpperLimit) l3OverFlow[3] ++;
	    break;
	  case 5 :
	    l3PtRate_motherBin_5->Fill(pt_L3);
	    if (pt_L3 > ptUpperLimit) l3OverFlow[4] ++;
	    break;
	  case 6 :
	    l3PtRate_motherBin_6->Fill(pt_L3);
	    if (pt_L3 > ptUpperLimit) l3OverFlow[5] ++;
	    break;
	  case 10 :
	    l3PtRate_motherBin_10->Fill(pt_L3);
	    if (pt_L3 > ptUpperLimit) l3OverFlow[9] ++;
	    break;
	  case 11 :
	    l3PtRate_motherBin_11->Fill(pt_L3);
	    if (pt_L3 > ptUpperLimit) l3OverFlow[10] ++;
	    break;
	  case 12 :
	    l3PtRate_motherBin_12->Fill(pt_L3);
	    if (pt_L3 > ptUpperLimit) l3OverFlow[11] ++;
	    break;
	  case 13 :
	    l3PtRate_motherBin_13->Fill(pt_L3);
	    if (pt_L3 > ptUpperLimit) l3OverFlow[12] ++;
	    break;
	  default :
	    //           cout << "oh, cock, this is interesting" << endl;
	    l3PtRate_motherBin_13->Fill(pt_L3);
	    if (pt_L3 > ptUpperLimit) l3OverFlow[12] ++;
	  }
	}
	else {
	  l3PtRate_motherBin_14->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3OverFlow[13] ++;
	}
      }
	else {
	  l3PtRate_motherBin_15->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3OverFlow[14] ++;
	}
      }
      
      if (passL2) {
	//double pt_L2 = findMaxPt(l2Pt,l2Eta,l2D0);
	//int l2_index = findIndexOfMaxPt(l2Pt,l2Eta,l2D0);
	int l2_index = iMu;
	double pt_L2 = (*l2Pt).at(iMu);
	if(!isData){
	if ((*l2IsAssociated).at(l2_index) == 1 && (*l2AssociationVar).at(l2_index) > -0.1) {
	  switch ( (*l2MotherBinNumber).at(l2_index) ) {
	  case 1 :
	    l2PtRate_motherBin_1->Fill(pt_L2);
	    if (pt_L2 > ptUpperLimit) l2OverFlow[0] ++;
	    break;
	  case 2 :
	    l2PtRate_motherBin_2->Fill(pt_L2);
	    if (pt_L2 > ptUpperLimit) l2OverFlow[1] ++;
	    break;
	  case 3 :
	    l2PtRate_motherBin_3->Fill(pt_L2);
	    if (pt_L2 > ptUpperLimit) l2OverFlow[2] ++;
	    break;
	  case 4 :
	    l2PtRate_motherBin_4->Fill(pt_L2);
	    if (pt_L2 > ptUpperLimit) l2OverFlow[3] ++;
	    break;
	  case 5 :
	    l2PtRate_motherBin_5->Fill(pt_L2);
	    if (pt_L2 > ptUpperLimit) l2OverFlow[4] ++;
	    break;
	  case 6 :
	    l2PtRate_motherBin_6->Fill(pt_L2);
	    if (pt_L2 > ptUpperLimit) l2OverFlow[5] ++;
	    break;
	  case 10 :
	    l2PtRate_motherBin_10->Fill(pt_L2);
	    if (pt_L2 > ptUpperLimit) l2OverFlow[9] ++;
	    break;
	  case 11 :
	    l2PtRate_motherBin_11->Fill(pt_L2);
	    if (pt_L2 > ptUpperLimit) l2OverFlow[10] ++;
	    break;
	  case 12 :
	    l2PtRate_motherBin_12->Fill(pt_L2);
	    if (pt_L2 > ptUpperLimit) l2OverFlow[11] ++;
	    break;
	  case 13 :
	    l2PtRate_motherBin_13->Fill(pt_L2);
	    if (pt_L2 > ptUpperLimit) l2OverFlow[12] ++;
	    break;
	  default :
	    //           cout << "oh, cock, this is interesting" << endl;
	    l2PtRate_motherBin_13->Fill(pt_L2);
	    if (pt_L2 > ptUpperLimit) l2OverFlow[12] ++;
	  }
	}
	else {
	  l2PtRate_motherBin_14->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2OverFlow[13] ++;
	}
	} else {
	  l2PtRate_motherBin_15->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2OverFlow[14] ++;
	}
      }
      
      if (passTK) {
	//double pt_tkTrack = findMaxPt(tkTrackPt,tkTrackEta,tkTrackD0);
	//int tkTrack_index = findIndexOfMaxPt(tkTrackPt,tkTrackEta,tkTrackD0);
	int tkTrack_index = iMu;
	double pt_tkTrack = (*tkTrackPt).at(iMu);
	if(!isData){
	if ((*tkTrackIsAssociated).at(tkTrack_index) == 1 && (*tkTrackAssociationVar).at(tkTrack_index) > -0.1 && fabs((*tkTrackAssociationPdgId).at(l3_index))==13) {
	  switch ( (*tkTrackMotherBinNumber).at(tkTrack_index) ) {
	  case 1 :
	    tkTrackPtRate_motherBin_1->Fill(pt_tkTrack);
	    if (pt_tkTrack > ptUpperLimit) tkTrackOverFlow[0] ++;
	    break;
	  case 2 :
	    tkTrackPtRate_motherBin_2->Fill(pt_tkTrack);
	    if (pt_tkTrack > ptUpperLimit) tkTrackOverFlow[1] ++;
	    break;
	  case 3 :
	    tkTrackPtRate_motherBin_3->Fill(pt_tkTrack);
	    if (pt_tkTrack > ptUpperLimit) tkTrackOverFlow[2] ++;
	    break;
	  case 4 :
	    tkTrackPtRate_motherBin_4->Fill(pt_tkTrack);
	    if (pt_tkTrack > ptUpperLimit) tkTrackOverFlow[3] ++;
	    break;
	  case 5 :
	    tkTrackPtRate_motherBin_5->Fill(pt_tkTrack);
	    if (pt_tkTrack > ptUpperLimit) tkTrackOverFlow[4] ++;
	    break;
	  case 6 :
	    tkTrackPtRate_motherBin_6->Fill(pt_tkTrack);
	    if (pt_tkTrack > ptUpperLimit) tkTrackOverFlow[5] ++;
	    break;
	  case 10 :
	    tkTrackPtRate_motherBin_10->Fill(pt_tkTrack);
	    if (pt_tkTrack > ptUpperLimit) tkTrackOverFlow[9] ++;
	    break;
	  case 11 :
	    tkTrackPtRate_motherBin_11->Fill(pt_tkTrack);
	    if (pt_tkTrack > ptUpperLimit) tkTrackOverFlow[10] ++;
	    break;
	  case 12 :
	    tkTrackPtRate_motherBin_12->Fill(pt_tkTrack);
	    if (pt_tkTrack > ptUpperLimit) tkTrackOverFlow[11] ++;
	    break;
	  case 13 :
	    tkTrackPtRate_motherBin_13->Fill(pt_tkTrack);
	    if (pt_tkTrack > ptUpperLimit) tkTrackOverFlow[12] ++;
	    break;
	  default :
	    //           cout << "oh, cock, this is interesting" << endl;
	    tkTrackPtRate_motherBin_13->Fill(pt_tkTrack);
	    if (pt_tkTrack > ptUpperLimit) tkTrackOverFlow[12] ++;
	  }
	}
	else {
	  tkTrackPtRate_motherBin_14->Fill(pt_tkTrack);
	  if (pt_tkTrack > ptUpperLimit) tkTrackOverFlow[13] ++;
	}
	} else {
	   tkTrackPtRate_motherBin_15->Fill(pt_tkTrack);
	   if (pt_tkTrack > ptUpperLimit) tkTrackOverFlow[14] ++; 
	}
      }
      
    }
  }
  
  fill_overflow(l3PtRate_motherBin_1,l3OverFlow[0]);
  fill_overflow(l3PtRate_motherBin_2,l3OverFlow[1]);
  fill_overflow(l3PtRate_motherBin_3,l3OverFlow[2]);
  fill_overflow(l3PtRate_motherBin_4,l3OverFlow[3]);
  fill_overflow(l3PtRate_motherBin_5,l3OverFlow[4]);
  fill_overflow(l3PtRate_motherBin_6,l3OverFlow[5]);
  fill_overflow(l3PtRate_motherBin_10,l3OverFlow[9]);
  fill_overflow(l3PtRate_motherBin_11,l3OverFlow[10]);
  fill_overflow(l3PtRate_motherBin_12,l3OverFlow[11]);
  fill_overflow(l3PtRate_motherBin_13,l3OverFlow[12]);
  fill_overflow(l3PtRate_motherBin_14,l3OverFlow[13]);
  fill_overflow(l3PtRate_motherBin_15,l3OverFlow[14]);

  make_cd_from_histo(l3PtRate_motherBin_1,l3PtRate_motherBin_1_cd);
  make_cd_from_histo(l3PtRate_motherBin_2,l3PtRate_motherBin_2_cd);
  make_cd_from_histo(l3PtRate_motherBin_3,l3PtRate_motherBin_3_cd);
  make_cd_from_histo(l3PtRate_motherBin_4,l3PtRate_motherBin_4_cd);
  make_cd_from_histo(l3PtRate_motherBin_5,l3PtRate_motherBin_5_cd);
  make_cd_from_histo(l3PtRate_motherBin_6,l3PtRate_motherBin_6_cd);
  make_cd_from_histo(l3PtRate_motherBin_10,l3PtRate_motherBin_10_cd);
  make_cd_from_histo(l3PtRate_motherBin_11,l3PtRate_motherBin_11_cd);
  make_cd_from_histo(l3PtRate_motherBin_12,l3PtRate_motherBin_12_cd);
  make_cd_from_histo(l3PtRate_motherBin_13,l3PtRate_motherBin_13_cd);
  make_cd_from_histo(l3PtRate_motherBin_14,l3PtRate_motherBin_14_cd);
  make_cd_from_histo(l3PtRate_motherBin_15,l3PtRate_motherBin_15_cd);

  if(!isData){
    l3PtRate->Add(l3PtRate_motherBin_1);
    l3PtRate->Add(l3PtRate_motherBin_2);
    //l3PtRate->Add(l3PtRate_motherBin_3);
    //l3PtRate->Add(l3PtRate_motherBin_4);
    //l3PtRate->Add(l3PtRate_motherBin_5);
    //l3PtRate->Add(l3PtRate_motherBin_6);
    l3PtRate->Add(l3PtRate_motherBin_10);
    l3PtRate->Add(l3PtRate_motherBin_11);
    l3PtRate->Add(l3PtRate_motherBin_12);
    l3PtRate->Add(l3PtRate_motherBin_13);
    l3PtRate->Add(l3PtRate_motherBin_14);
  } else {
    l3PtRate->Add(l3PtRate_motherBin_15);
  }

  if(!isData){
    l3PtRate_cd->Add(l3PtRate_motherBin_1_cd);
    l3PtRate_cd->Add(l3PtRate_motherBin_2_cd);
    //l3PtRate_cd->Add(l3PtRate_motherBin_3_cd);
    //l3PtRate_cd->Add(l3PtRate_motherBin_4_cd);
    //l3PtRate_cd->Add(l3PtRate_motherBin_5_cd);
    //l3PtRate_cd->Add(l3PtRate_motherBin_6_cd);
    l3PtRate_cd->Add(l3PtRate_motherBin_10_cd);
    l3PtRate_cd->Add(l3PtRate_motherBin_11_cd);
    l3PtRate_cd->Add(l3PtRate_motherBin_12_cd);
    l3PtRate_cd->Add(l3PtRate_motherBin_13_cd);
    l3PtRate_cd->Add(l3PtRate_motherBin_14_cd);
  } else {
    l3PtRate_cd->Add(l3PtRate_motherBin_15_cd);
  }
  //aaa
  
  fill_overflow(l2PtRate_motherBin_1,l2OverFlow[0]);
  fill_overflow(l2PtRate_motherBin_2,l2OverFlow[1]);
  fill_overflow(l2PtRate_motherBin_3,l2OverFlow[2]);
  fill_overflow(l2PtRate_motherBin_4,l2OverFlow[3]);
  fill_overflow(l2PtRate_motherBin_5,l2OverFlow[4]);
  fill_overflow(l2PtRate_motherBin_6,l2OverFlow[5]);
  fill_overflow(l2PtRate_motherBin_10,l2OverFlow[9]);
  fill_overflow(l2PtRate_motherBin_11,l2OverFlow[10]);
  fill_overflow(l2PtRate_motherBin_12,l2OverFlow[11]);
  fill_overflow(l2PtRate_motherBin_13,l2OverFlow[12]);
  fill_overflow(l2PtRate_motherBin_14,l2OverFlow[13]);
  fill_overflow(l2PtRate_motherBin_15,l2OverFlow[14]);

  make_cd_from_histo(l2PtRate_motherBin_1,l2PtRate_motherBin_1_cd);
  make_cd_from_histo(l2PtRate_motherBin_2,l2PtRate_motherBin_2_cd);
  make_cd_from_histo(l2PtRate_motherBin_3,l2PtRate_motherBin_3_cd);
  make_cd_from_histo(l2PtRate_motherBin_4,l2PtRate_motherBin_4_cd);
  make_cd_from_histo(l2PtRate_motherBin_5,l2PtRate_motherBin_5_cd);
  make_cd_from_histo(l2PtRate_motherBin_6,l2PtRate_motherBin_6_cd);
  make_cd_from_histo(l2PtRate_motherBin_10,l2PtRate_motherBin_10_cd);
  make_cd_from_histo(l2PtRate_motherBin_11,l2PtRate_motherBin_11_cd);
  make_cd_from_histo(l2PtRate_motherBin_12,l2PtRate_motherBin_12_cd);
  make_cd_from_histo(l2PtRate_motherBin_13,l2PtRate_motherBin_13_cd);
  make_cd_from_histo(l2PtRate_motherBin_14,l2PtRate_motherBin_14_cd);
  make_cd_from_histo(l2PtRate_motherBin_15,l2PtRate_motherBin_15_cd);

  if(!isData){
  l2PtRate->Add(l2PtRate_motherBin_1);
  l2PtRate->Add(l2PtRate_motherBin_2);
  //l2PtRate->Add(l2PtRate_motherBin_3);
  //l2PtRate->Add(l2PtRate_motherBin_4);
  //l2PtRate->Add(l2PtRate_motherBin_5);
  //l2PtRate->Add(l2PtRate_motherBin_6);
  l2PtRate->Add(l2PtRate_motherBin_10);
  l2PtRate->Add(l2PtRate_motherBin_11);
  l2PtRate->Add(l2PtRate_motherBin_12);
  l2PtRate->Add(l2PtRate_motherBin_13);
  l2PtRate->Add(l2PtRate_motherBin_14);
  } else {
    l2PtRate->Add(l2PtRate_motherBin_15);
  }

  if(!isData){
    l2PtRate_cd->Add(l2PtRate_motherBin_1_cd);
    l2PtRate_cd->Add(l2PtRate_motherBin_2_cd);
    //l2PtRate_cd->Add(l2PtRate_motherBin_3_cd);
    //l2PtRate_cd->Add(l2PtRate_motherBin_4_cd);
    //l2PtRate_cd->Add(l2PtRate_motherBin_5_cd);
    //l2PtRate_cd->Add(l2PtRate_motherBin_6_cd);
    l2PtRate_cd->Add(l2PtRate_motherBin_10_cd);
    l2PtRate_cd->Add(l2PtRate_motherBin_11_cd);
    l2PtRate_cd->Add(l2PtRate_motherBin_12_cd);
    l2PtRate_cd->Add(l2PtRate_motherBin_13_cd);
    l2PtRate_cd->Add(l2PtRate_motherBin_14_cd);
  } else {
    l2PtRate_cd->Add(l2PtRate_motherBin_15_cd);
  }
  //aaa

  fill_overflow(tkTrackPtRate_motherBin_1,tkTrackOverFlow[0]);
  fill_overflow(tkTrackPtRate_motherBin_2,tkTrackOverFlow[1]);
  fill_overflow(tkTrackPtRate_motherBin_3,tkTrackOverFlow[2]);
  fill_overflow(tkTrackPtRate_motherBin_4,tkTrackOverFlow[3]);
  fill_overflow(tkTrackPtRate_motherBin_5,tkTrackOverFlow[4]);
  fill_overflow(tkTrackPtRate_motherBin_6,tkTrackOverFlow[5]);
  fill_overflow(tkTrackPtRate_motherBin_10,tkTrackOverFlow[9]);
  fill_overflow(tkTrackPtRate_motherBin_11,tkTrackOverFlow[10]);
  fill_overflow(tkTrackPtRate_motherBin_12,tkTrackOverFlow[11]);
  fill_overflow(tkTrackPtRate_motherBin_13,tkTrackOverFlow[12]);
  fill_overflow(tkTrackPtRate_motherBin_14,tkTrackOverFlow[13]);
  fill_overflow(tkTrackPtRate_motherBin_15,tkTrackOverFlow[14]);

  make_cd_from_histo(tkTrackPtRate_motherBin_1,tkTrackPtRate_motherBin_1_cd);
  make_cd_from_histo(tkTrackPtRate_motherBin_2,tkTrackPtRate_motherBin_2_cd);
  make_cd_from_histo(tkTrackPtRate_motherBin_3,tkTrackPtRate_motherBin_3_cd);
  make_cd_from_histo(tkTrackPtRate_motherBin_4,tkTrackPtRate_motherBin_4_cd);
  make_cd_from_histo(tkTrackPtRate_motherBin_5,tkTrackPtRate_motherBin_5_cd);
  make_cd_from_histo(tkTrackPtRate_motherBin_6,tkTrackPtRate_motherBin_6_cd);
  make_cd_from_histo(tkTrackPtRate_motherBin_10,tkTrackPtRate_motherBin_10_cd);
  make_cd_from_histo(tkTrackPtRate_motherBin_11,tkTrackPtRate_motherBin_11_cd);
  make_cd_from_histo(tkTrackPtRate_motherBin_12,tkTrackPtRate_motherBin_12_cd);
  make_cd_from_histo(tkTrackPtRate_motherBin_13,tkTrackPtRate_motherBin_13_cd);
  make_cd_from_histo(tkTrackPtRate_motherBin_14,tkTrackPtRate_motherBin_14_cd);
  make_cd_from_histo(tkTrackPtRate_motherBin_15,tkTrackPtRate_motherBin_15_cd);

  if(!isData){
    tkTrackPtRate->Add(tkTrackPtRate_motherBin_1);
    tkTrackPtRate->Add(tkTrackPtRate_motherBin_2);
    //tkTrackPtRate->Add(tkTrackPtRate_motherBin_3);
    //tkTrackPtRate->Add(tkTrackPtRate_motherBin_4);
    //tkTrackPtRate->Add(tkTrackPtRate_motherBin_5);
    //tkTrackPtRate->Add(tkTrackPtRate_motherBin_6);
    tkTrackPtRate->Add(tkTrackPtRate_motherBin_10);
    tkTrackPtRate->Add(tkTrackPtRate_motherBin_11);
    tkTrackPtRate->Add(tkTrackPtRate_motherBin_12);
    tkTrackPtRate->Add(tkTrackPtRate_motherBin_13);
    tkTrackPtRate->Add(tkTrackPtRate_motherBin_14);
  } else {
    tkTrackPtRate->Add(tkTrackPtRate_motherBin_15);
  }
  
  if(!isData){
    tkTrackPtRate_cd->Add(tkTrackPtRate_motherBin_1_cd);
    tkTrackPtRate_cd->Add(tkTrackPtRate_motherBin_2_cd);
    //tkTrackPtRate_cd->Add(tkTrackPtRate_motherBin_3_cd);
    //tkTrackPtRate_cd->Add(tkTrackPtRate_motherBin_4_cd);
    //tkTrackPtRate_cd->Add(tkTrackPtRate_motherBin_5_cd);
    //tkTrackPtRate_cd->Add(tkTrackPtRate_motherBin_6_cd);
    tkTrackPtRate_cd->Add(tkTrackPtRate_motherBin_10_cd);
    tkTrackPtRate_cd->Add(tkTrackPtRate_motherBin_11_cd);
    tkTrackPtRate_cd->Add(tkTrackPtRate_motherBin_12_cd);
    tkTrackPtRate_cd->Add(tkTrackPtRate_motherBin_13_cd);
    tkTrackPtRate_cd->Add(tkTrackPtRate_motherBin_14_cd);
  } else {
    tkTrackPtRate_cd->Add(tkTrackPtRate_motherBin_15_cd);
  }
  //aaaa

  //  double crossSection_mb = 75.28; // MinBias
  double crossSection_mb = 51.6; // ppMuX
  //  double crossSection_mb = 0.000317 ; // ttbar
  //  double filterEff = 1.0; //MinBias
  double filterEff = 0.00289; // ppMuX
  //   double filterEff = 0.33;  // ttbar

  double eventNumberUnit= 1000000;
  //   double filterEff = inclusiveppMuX_filterEff;
  double numberOfEvents = nEntries /eventNumberUnit /filterEff ;
  cout << "numberOfEvents = " << numberOfEvents << endl;

  double crossSection=crossSection_mb ;
  double L_E32= 100000; //10^32 cm-2.s-1 = 10^5 mb-1.s-1
  double L_E31= 10000; //10^31 cm-2.s-1 = 10^4 mb-1.s-1
  double L_E30= 1000; //10^30 cm-2.s-1 = 10^3 mb-1.s-1
  double L=L_E30;

  double mbInvToHz=L/eventNumberUnit;
  //  double rateFactor=crossSection / numberOfEvents * mbInvToHz ;
  double rateFactorMC = 0.048492353331627547;
  rateFactorMC = rateFactorMC * 4;

  double rateFactor = (weight==0) ? rateFactorMC : weight;

  cout << "rateFactor = " << rateFactor << endl;

  histDir->cd();
  scaleToRate("*PtRate*", rateFactor);


  histDir->Write("",TObject::kOverwrite);

  histFile->Close();
  
  return 0;
}




double findMaxPt(vector<double> *pts) { // just a simple "find max"
  double maxPt = -999;
  for (int i = 0; i < pts->size(); i++) {
    if (pts->at(i) > maxPt) maxPt = pts->at(i);
  }
  return maxPt;
}

double findMaxPt(vector<double> *pts, vector<double> *etas) { // find max for l2
  double maxPt = -999;
  double etaCut = 2.5;
  for (int i = 0; i < pts->size(); i++) {
    if (pts->at(i) > maxPt && fabs(etas->at(i)) < etaCut) maxPt = pts->at(i);
  }
  return maxPt;
}

double findMaxPt(vector<double> *pts, vector<double> *etas, vector<double> *d0s) { // find max for l3
  double maxPt = -999;
  double etaCut = 1000000000; // junk, just to make sure I get no fucking seg faults
  double d0Cut = 2;
  for (int i = 0; i < pts->size(); i++) {
    if (pts->at(i) > maxPt && fabs(etas->at(i)) < etaCut && fabs(d0s->at(i)) < d0Cut) maxPt = pts->at(i);
  }
  return maxPt;
}

int findIndexOfMaxPt(vector<double> *pts) { // just a simple "find max"
  double maxPt = -999;
  int index = -999;
  for (int i = 0; i < pts->size(); i++) {
    if (pts->at(i) > maxPt) {
      index = i;
      maxPt = pts->at(i);
    }
  }
  return index;
}

int findIndexOfMaxPt(vector<double> *pts, vector<double> *etas) { // find max for l2
  double maxPt = -999;
  double etaCut = 2.5;
  int index = -999;
  for (int i = 0; i < pts->size(); i++) {
    if (pts->at(i) > maxPt && fabs(etas->at(i)) < etaCut) {
      index = i;
      maxPt = pts->at(i);
    }
  }
  return index;
}

int findIndexOfMaxPt(vector<double> *pts, vector<double> *etas, vector<double> *d0s) { // find max for l3
  double maxPt = -999;
  double etaCut = 999;
  double d0Cut = 2;
  int index = -999;
  for (int i = 0; i < pts->size(); i++) {
    if (pts->at(i) > maxPt && fabs(etas->at(i)) < etaCut && fabs(d0s->at(i)) < d0Cut) {
      index = i;
      maxPt = pts->at(i);
    }
  }
  return index;
}

void fill_overflow(TH1F* hist, double overflow) {
  int n_max = hist->GetNbinsX();
  hist->SetBinContent(n_max + 1, overflow);
}

// This one takes a hist and makes a cd.  Doesn't have to be same number of bins, but the number of bins should be the same
void make_cd_from_histo(TH1F* hist, TH1F* cd) {

  int nBins1 = hist->GetNbinsX();
  int nBins2 = cd->GetNbinsX();
  if (hist->GetBinWidth(1) != cd->GetBinWidth(1)) {
    cout << "This is wrong." << endl;
    cout << "width 1, 2 = " << hist->GetBinWidth(1) << " " << cd->GetBinWidth(1) << endl;
    return;
  }

  double content = 0;
  for(int i=1; i<=nBins2; ++i){
    for (int j=i; j<=nBins1 + 1; ++j) { // include overflow
      content+=hist->GetBinContent(j);
    }
    cd->SetBinContent(i, content); //Integrate.
    content = 0; // Reset it, you muppet.
  }

  /*  double norm_factor = hist->GetEntries();
      
  if (norm_factor != 0) {
  cd->Scale(1/norm_factor);
  hist->Scale(1/norm_factor);
  }
  */
  return;
}

void scaleToRate(const char* collectionName, double rateFactor) {
  TRegexp reg(collectionName, kTRUE);
  
  //    gDirectory->ls();
  TList* list = gDirectory->GetList() ;
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  
  while (obj = iter->Next()) {
    if (! obj->InheritsFrom(TH1::Class())) {
      //      cout << "bugger" << endl;
      continue;
    }
    
    
    TString name = obj->GetName();
    cout << "Testing name: " << name << " against " << collectionName << endl;
    
    if (TString(collectionName).MaybeRegexp()) {
      cout << "we have a possible match" << endl;
      cout << "Trying to match to " << TString(obj->GetName()) << endl;
      if (TString(obj->GetName()).Index(reg) < 0 ) {
	cout << "failure here.  Argument returns " << TString(obj->GetName()).Index(reg) << endl;
	continue;
      }
    }
    else if (! name.BeginsWith(collectionName)) continue;
    
    cout << "We're trying to scale" << name << endl;
    ((TH1*)obj)->Scale(rateFactor);
    
  }
  
}
