#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
//#include "Math/LorentzVector.h"
#include <vector>
#include "RecoTree.h"

//int ScanTree ( TTree* tree) {
int ScanTreeD0 ( TTree* tree, char *fileName, bool isData=false, double weight = 1.0) {

  // This reads in the tree.  As you might imagine.
  Init(tree);
  TFile *histFile = new TFile(fileName,"UPDATE");
  TDirectory *histDir = histFile->mkdir("D0Rate");

  histDir->cd();

  // Time for some histograms and stacks

  double lowerLimit = -5;
  double upperLimit = 5;
  double int nBins    = 50;

  TH1F *l3D0Rate_motherBin_1 = new TH1F("l3D0Rate_motherBin_1"," #pi^{+/-}",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_2 = new TH1F("l3D0Rate_motherBin_2"," K",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_3 = new TH1F("l3D0Rate_motherBin_3"," D",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_4 = new TH1F("l3D0Rate_motherBin_4"," B",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_5 = new TH1F("l3D0Rate_motherBin_5"," #Lambda_{b}",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_6 = new TH1F("l3D0Rate_motherBin_6"," J/#Psi",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_10 = new TH1F("l3D0Rate_motherBin_10"," #tau^{+/-}",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_11 = new TH1F("l3D0Rate_motherBin_11"," #mu",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_12 = new TH1F("l3D0Rate_motherBin_12"," b/c",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_13 = new TH1F("l3D0Rate_motherBin_13"," other",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_14 = new TH1F("l3D0Rate_motherBin_14"," non-associated",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_15 = new TH1F("l3D0Rate_motherBin_15"," Data",nBins,lowerLimit,upperLimit);

  TH1F *l3D0Rate_motherBin_1_cd = new TH1F("l3D0Rate_motherBin_1_cd"," #pi^{+/-}",nBins,lowerLimit,upperLimit);  
  TH1F *l3D0Rate_motherBin_2_cd = new TH1F("l3D0Rate_motherBin_2_cd"," K",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_3_cd = new TH1F("l3D0Rate_motherBin_3_cd"," D",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_4_cd = new TH1F("l3D0Rate_motherBin_4_cd"," B",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_5_cd = new TH1F("l3D0Rate_motherBin_5_cd"," #Lambda_{b}",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_6_cd = new TH1F("l3D0Rate_motherBin_6_cd"," J/#Psi",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_10_cd = new TH1F("l3D0Rate_motherBin_10_cd"," #tau^{+/-}",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_11_cd = new TH1F("l3D0Rate_motherBin_11_cd"," #mu",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_12_cd = new TH1F("l3D0Rate_motherBin_12_cd"," b/c",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_13_cd = new TH1F("l3D0Rate_motherBin_13_cd"," other",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_14_cd = new TH1F("l3D0Rate_motherBin_14_cd"," non-associated",nBins,lowerLimit,upperLimit);
  TH1F *l3D0Rate_motherBin_15_cd = new TH1F("l3D0Rate_motherBin_15_cd"," Data",nBins,lowerLimit,upperLimit);



  THStack *l3D0Rate = new THStack("l3D0Rate","L3 rate as f(p_{T,L3})");
  THStack *l3D0Rate_cd = new THStack("l3D0Rate_cd","L3 rate as f(p_{T,L3})");
  //  THStack *l3PreIsoD0Rate = new THStack("l3PreIsoD0Rate","L3PreIso rate as f(p_{T,L3})");
  //  THStack *l3PreIsoD0Rate_cd = new THStack("l3PreIsoD0Rate_cd","L3PreIso rate as f(p_{T,L3})");
  //  THStack *l3IsoD0Rate = new THStack("l3IsoD0Rate","L3Iso rate as f(p_{T,L3})");
  //  THStack *l3IsoD0Rate_cd = new THStack("l3IsoD0Rate_cd","L3Iso rate as f(p_{T,L3})");

  TH1F *l2D0Rate_motherBin_1 = new TH1F("l2D0Rate_motherBin_1"," #pi^{+/-}",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_2 = new TH1F("l2D0Rate_motherBin_2"," K",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_3 = new TH1F("l2D0Rate_motherBin_3"," D",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_4 = new TH1F("l2D0Rate_motherBin_4"," B",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_5 = new TH1F("l2D0Rate_motherBin_5"," #Lambda_{b}",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_6 = new TH1F("l2D0Rate_motherBin_6"," J/#Psi",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_10 = new TH1F("l2D0Rate_motherBin_10"," #tau^{+/-}",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_11 = new TH1F("l2D0Rate_motherBin_11"," #mu",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_12 = new TH1F("l2D0Rate_motherBin_12"," b/c",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_13 = new TH1F("l2D0Rate_motherBin_13"," other",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_14 = new TH1F("l2D0Rate_motherBin_14"," non-associated",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_15 = new TH1F("l2D0Rate_motherBin_15"," Data",nBins,lowerLimit,upperLimit);

  TH1F *l2D0Rate_motherBin_1_cd = new TH1F("l2D0Rate_motherBin_1_cd"," #pi^{+/-}",nBins,lowerLimit,upperLimit);  
  TH1F *l2D0Rate_motherBin_2_cd = new TH1F("l2D0Rate_motherBin_2_cd"," K",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_3_cd = new TH1F("l2D0Rate_motherBin_3_cd"," D",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_4_cd = new TH1F("l2D0Rate_motherBin_4_cd"," B",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_5_cd = new TH1F("l2D0Rate_motherBin_5_cd"," #Lambda_{b}",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_6_cd = new TH1F("l2D0Rate_motherBin_6_cd"," J/#Psi",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_10_cd = new TH1F("l2D0Rate_motherBin_10_cd"," #tau^{+/-}",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_11_cd = new TH1F("l2D0Rate_motherBin_11_cd"," #mu",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_12_cd = new TH1F("l2D0Rate_motherBin_12_cd"," b/c",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_13_cd = new TH1F("l2D0Rate_motherBin_13_cd"," other",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_14_cd = new TH1F("l2D0Rate_motherBin_14_cd"," non-associated",nBins,lowerLimit,upperLimit);
  TH1F *l2D0Rate_motherBin_15_cd = new TH1F("l2D0Rate_motherBin_15_cd"," Data",nBins,lowerLimit,upperLimit);



  THStack *l2D0Rate = new THStack("l2D0Rate","L2 rate as f(p_{T,L2})");
  THStack *l2D0Rate_cd = new THStack("l2D0Rate_cd","L2 rate as f(p_{T,L2})");
  //  THStack *l2PreIsoD0Rate = new THStack("l2PreIsoD0Rate","L2PreIso rate as f(p_{T,L2})");
  //  THStack *l2PreIsoD0Rate_cd = new THStack("l2PreIsoD0Rate_cd","L2PreIso rate as f(p_{T,L2})");
  //  THStack *l2IsoD0Rate = new THStack("l2IsoD0Rate","L2Iso rate as f(p_{T,L2})");
  //  THStack *l2IsoD0Rate_cd = new THStack("l2IsoD0Rate_cd","L2Iso rate as f(p_{T,L2})");

  TH1F *tkTrackD0Rate_motherBin_1 = new TH1F("tkTrackD0Rate_motherBin_1"," #pi^{+/-}",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_2 = new TH1F("tkTrackD0Rate_motherBin_2"," K",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_3 = new TH1F("tkTrackD0Rate_motherBin_3"," D",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_4 = new TH1F("tkTrackD0Rate_motherBin_4"," B",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_5 = new TH1F("tkTrackD0Rate_motherBin_5"," #Lambda_{b}",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_6 = new TH1F("tkTrackD0Rate_motherBin_6"," J/#Psi",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_10 = new TH1F("tkTrackD0Rate_motherBin_10"," #tau^{+/-}",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_11 = new TH1F("tkTrackD0Rate_motherBin_11"," #mu",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_12 = new TH1F("tkTrackD0Rate_motherBin_12"," b/c",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_13 = new TH1F("tkTrackD0Rate_motherBin_13"," other",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_14 = new TH1F("tkTrackD0Rate_motherBin_14"," non-associated",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_15 = new TH1F("tkTrackD0Rate_motherBin_15"," Data",nBins,lowerLimit,upperLimit);

  TH1F *tkTrackD0Rate_motherBin_1_cd = new TH1F("tkTrackD0Rate_motherBin_1_cd"," #pi^{+/-}",nBins,lowerLimit,upperLimit);  
  TH1F *tkTrackD0Rate_motherBin_2_cd = new TH1F("tkTrackD0Rate_motherBin_2_cd"," K",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_3_cd = new TH1F("tkTrackD0Rate_motherBin_3_cd"," D",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_4_cd = new TH1F("tkTrackD0Rate_motherBin_4_cd"," B",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_5_cd = new TH1F("tkTrackD0Rate_motherBin_5_cd"," #Lambda_{b}",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_6_cd = new TH1F("tkTrackD0Rate_motherBin_6_cd"," J/#Psi",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_10_cd = new TH1F("tkTrackD0Rate_motherBin_10_cd"," #tau^{+/-}",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_11_cd = new TH1F("tkTrackD0Rate_motherBin_11_cd"," #mu",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_12_cd = new TH1F("tkTrackD0Rate_motherBin_12_cd"," b/c",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_13_cd = new TH1F("tkTrackD0Rate_motherBin_13_cd"," other",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_14_cd = new TH1F("tkTrackD0Rate_motherBin_14_cd"," non-associated",nBins,lowerLimit,upperLimit);
  TH1F *tkTrackD0Rate_motherBin_15_cd = new TH1F("tkTrackD0Rate_motherBin_15_cd"," Data",nBins,lowerLimit,upperLimit);



  THStack *tkTrackD0Rate = new THStack("tkTrackD0Rate","tkTrack rate as f(p_{T,tkTrack})");
  THStack *tkTrackD0Rate_cd = new THStack("tkTrackD0Rate_cd","tkTrack rate as f(p_{T,tkTrack})");
  //  THStack *tkTrackPreIsoD0Rate = new THStack("tkTrackPreIsoD0Rate","tkTrackPreIso rate as f(p_{T,tkTrack})");
  //  THStack *tkTrackPreIsoD0Rate_cd = new THStack("tkTrackPreIsoD0Rate_cd","tkTrackPreIso rate as f(p_{T,tkTrack})");
  //  THStack *tkTrackIsoD0Rate = new THStack("tkTrackIsoD0Rate","tkTrackIso rate as f(p_{T,tkTrack})");
  //  THStack *tkTrackIsoD0Rate_cd = new THStack("tkTrackIsoD0Rate_cd","tkTrackIso rate as f(p_{T,tkTrack})");

  // There are many comments to be made about ROOT.  Most of them Rated-R or Rated-X.
  // Here, I'll limit myself to saying this is necessary for the TStacks to go into the TDirectory.
  gDirectory->Append(l2D0Rate);
  gDirectory->Append(l3D0Rate);
  gDirectory->Append(tkTrackD0Rate);

  gDirectory->Append(l2D0Rate_cd);
  gDirectory->Append(l3D0Rate_cd);
  gDirectory->Append(tkTrackD0Rate_cd);

  //  gDirectory->Append(l2IsoD0Rate);
  //  gDirectory->Append(l3PreIsoD0Rate);
  //  gDirectory->Append(l3IsoD0Rate);

  //  gDirectory->Append(l2IsoD0Rate_cd);
  //  gDirectory->Append(l3PreIsoD0Rate_cd);
  //  gDirectory->Append(l3IsoD0Rate_cd);

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
	double pt_L3 = (*l3D0).at(iMu);
	if(!isData){
	if ((*l3IsAssociated).at(l3_index) == 1 && (*l3AssociationVar).at(l3_index) > -0.1  && fabs((*l3AssociationPdgId).at(l3_index))==13) {
	  switch ( (*l3MotherBinNumber).at(l3_index) ) {
	  case 1 :
	    l3D0Rate_motherBin_1->Fill(pt_L3);
	    if (pt_L3 > upperLimit) l3OverFlow[0] ++;
	    break;
	  case 2 :
	    l3D0Rate_motherBin_2->Fill(pt_L3);
	    if (pt_L3 > upperLimit) l3OverFlow[1] ++;
	    break;
	  case 3 :
	    l3D0Rate_motherBin_3->Fill(pt_L3);
	    if (pt_L3 > upperLimit) l3OverFlow[2] ++;
	    break;
	  case 4 :
	    l3D0Rate_motherBin_4->Fill(pt_L3);
	    if (pt_L3 > upperLimit) l3OverFlow[3] ++;
	    break;
	  case 5 :
	    l3D0Rate_motherBin_5->Fill(pt_L3);
	    if (pt_L3 > upperLimit) l3OverFlow[4] ++;
	    break;
	  case 6 :
	    l3D0Rate_motherBin_6->Fill(pt_L3);
	    if (pt_L3 > upperLimit) l3OverFlow[5] ++;
	    break;
	  case 10 :
	    l3D0Rate_motherBin_10->Fill(pt_L3);
	    if (pt_L3 > upperLimit) l3OverFlow[9] ++;
	    break;
	  case 11 :
	    l3D0Rate_motherBin_11->Fill(pt_L3);
	    if (pt_L3 > upperLimit) l3OverFlow[10] ++;
	    break;
	  case 12 :
	    l3D0Rate_motherBin_12->Fill(pt_L3);
	    if (pt_L3 > upperLimit) l3OverFlow[11] ++;
	    break;
	  case 13 :
	    l3D0Rate_motherBin_13->Fill(pt_L3);
	    if (pt_L3 > upperLimit) l3OverFlow[12] ++;
	    break;
	  default :
	    //           cout << "oh, cock, this is interesting" << endl;
	    l3D0Rate_motherBin_13->Fill(pt_L3);
	    if (pt_L3 > upperLimit) l3OverFlow[12] ++;
	  }
	}
	else {
	  l3D0Rate_motherBin_14->Fill(pt_L3);
	  if (pt_L3 > upperLimit) l3OverFlow[13] ++;
	}
      }
	else {
	  l3D0Rate_motherBin_15->Fill(pt_L3);
	  if (pt_L3 > upperLimit) l3OverFlow[14] ++;
	}
      }
      
      if (passL2) {
	//double pt_L2 = findMaxPt(l2Pt,l2Eta,l2D0);
	//int l2_index = findIndexOfMaxPt(l2Pt,l2Eta,l2D0);
	int l2_index = iMu;
	double pt_L2 = (*l2D0).at(iMu);
	if(!isData){
	if ((*l2IsAssociated).at(l2_index) == 1 && (*l2AssociationVar).at(l2_index) > -0.1) {
	  switch ( (*l2MotherBinNumber).at(l2_index) ) {
	  case 1 :
	    l2D0Rate_motherBin_1->Fill(pt_L2);
	    if (pt_L2 > upperLimit) l2OverFlow[0] ++;
	    break;
	  case 2 :
	    l2D0Rate_motherBin_2->Fill(pt_L2);
	    if (pt_L2 > upperLimit) l2OverFlow[1] ++;
	    break;
	  case 3 :
	    l2D0Rate_motherBin_3->Fill(pt_L2);
	    if (pt_L2 > upperLimit) l2OverFlow[2] ++;
	    break;
	  case 4 :
	    l2D0Rate_motherBin_4->Fill(pt_L2);
	    if (pt_L2 > upperLimit) l2OverFlow[3] ++;
	    break;
	  case 5 :
	    l2D0Rate_motherBin_5->Fill(pt_L2);
	    if (pt_L2 > upperLimit) l2OverFlow[4] ++;
	    break;
	  case 6 :
	    l2D0Rate_motherBin_6->Fill(pt_L2);
	    if (pt_L2 > upperLimit) l2OverFlow[5] ++;
	    break;
	  case 10 :
	    l2D0Rate_motherBin_10->Fill(pt_L2);
	    if (pt_L2 > upperLimit) l2OverFlow[9] ++;
	    break;
	  case 11 :
	    l2D0Rate_motherBin_11->Fill(pt_L2);
	    if (pt_L2 > upperLimit) l2OverFlow[10] ++;
	    break;
	  case 12 :
	    l2D0Rate_motherBin_12->Fill(pt_L2);
	    if (pt_L2 > upperLimit) l2OverFlow[11] ++;
	    break;
	  case 13 :
	    l2D0Rate_motherBin_13->Fill(pt_L2);
	    if (pt_L2 > upperLimit) l2OverFlow[12] ++;
	    break;
	  default :
	    //           cout << "oh, cock, this is interesting" << endl;
	    l2D0Rate_motherBin_13->Fill(pt_L2);
	    if (pt_L2 > upperLimit) l2OverFlow[12] ++;
	  }
	}
	else {
	  l2D0Rate_motherBin_14->Fill(pt_L2);
	  if (pt_L2 > upperLimit) l2OverFlow[13] ++;
	}
	} else {
	  l2D0Rate_motherBin_15->Fill(pt_L2);
	  if (pt_L2 > upperLimit) l2OverFlow[14] ++;
	}
      }
      
      if (passTK) {
	//double pt_tkTrack = findMaxPt(tkTrackPt,tkTrackEta,tkTrackD0);
	//int tkTrack_index = findIndexOfMaxPt(tkTrackPt,tkTrackEta,tkTrackD0);
	int tkTrack_index = iMu;
	double pt_tkTrack = (*tkTrackD0).at(iMu);
	if(!isData){
	if ((*tkTrackIsAssociated).at(tkTrack_index) == 1 && (*tkTrackAssociationVar).at(tkTrack_index) > -0.1  && fabs((*tkTrackAssociationPdgId).at(l3_index))==13) {
	  switch ( (*tkTrackMotherBinNumber).at(tkTrack_index) ) {
	  case 1 :
	    tkTrackD0Rate_motherBin_1->Fill(pt_tkTrack);
	    if (pt_tkTrack > upperLimit) tkTrackOverFlow[0] ++;
	    break;
	  case 2 :
	    tkTrackD0Rate_motherBin_2->Fill(pt_tkTrack);
	    if (pt_tkTrack > upperLimit) tkTrackOverFlow[1] ++;
	    break;
	  case 3 :
	    tkTrackD0Rate_motherBin_3->Fill(pt_tkTrack);
	    if (pt_tkTrack > upperLimit) tkTrackOverFlow[2] ++;
	    break;
	  case 4 :
	    tkTrackD0Rate_motherBin_4->Fill(pt_tkTrack);
	    if (pt_tkTrack > upperLimit) tkTrackOverFlow[3] ++;
	    break;
	  case 5 :
	    tkTrackD0Rate_motherBin_5->Fill(pt_tkTrack);
	    if (pt_tkTrack > upperLimit) tkTrackOverFlow[4] ++;
	    break;
	  case 6 :
	    tkTrackD0Rate_motherBin_6->Fill(pt_tkTrack);
	    if (pt_tkTrack > upperLimit) tkTrackOverFlow[5] ++;
	    break;
	  case 10 :
	    tkTrackD0Rate_motherBin_10->Fill(pt_tkTrack);
	    if (pt_tkTrack > upperLimit) tkTrackOverFlow[9] ++;
	    break;
	  case 11 :
	    tkTrackD0Rate_motherBin_11->Fill(pt_tkTrack);
	    if (pt_tkTrack > upperLimit) tkTrackOverFlow[10] ++;
	    break;
	  case 12 :
	    tkTrackD0Rate_motherBin_12->Fill(pt_tkTrack);
	    if (pt_tkTrack > upperLimit) tkTrackOverFlow[11] ++;
	    break;
	  case 13 :
	    tkTrackD0Rate_motherBin_13->Fill(pt_tkTrack);
	    if (pt_tkTrack > upperLimit) tkTrackOverFlow[12] ++;
	    break;
	  default :
	    //           cout << "oh, cock, this is interesting" << endl;
	    tkTrackD0Rate_motherBin_13->Fill(pt_tkTrack);
	    if (pt_tkTrack > upperLimit) tkTrackOverFlow[12] ++;
	  }
	}
	else {
	  tkTrackD0Rate_motherBin_14->Fill(pt_tkTrack);
	  if (pt_tkTrack > upperLimit) tkTrackOverFlow[13] ++;
	}
	} else {
	   tkTrackD0Rate_motherBin_15->Fill(pt_tkTrack);
	   if (pt_tkTrack > upperLimit) tkTrackOverFlow[14] ++; 
	}
      }
      
    }
  }
  
  fill_overflow(l3D0Rate_motherBin_1,l3OverFlow[0]);
  fill_overflow(l3D0Rate_motherBin_2,l3OverFlow[1]);
  fill_overflow(l3D0Rate_motherBin_3,l3OverFlow[2]);
  fill_overflow(l3D0Rate_motherBin_4,l3OverFlow[3]);
  fill_overflow(l3D0Rate_motherBin_5,l3OverFlow[4]);
  fill_overflow(l3D0Rate_motherBin_6,l3OverFlow[5]);
  fill_overflow(l3D0Rate_motherBin_10,l3OverFlow[9]);
  fill_overflow(l3D0Rate_motherBin_11,l3OverFlow[10]);
  fill_overflow(l3D0Rate_motherBin_12,l3OverFlow[11]);
  fill_overflow(l3D0Rate_motherBin_13,l3OverFlow[12]);
  fill_overflow(l3D0Rate_motherBin_14,l3OverFlow[13]);
  fill_overflow(l3D0Rate_motherBin_15,l3OverFlow[14]);

  make_cd_from_histo(l3D0Rate_motherBin_1,l3D0Rate_motherBin_1_cd);
  make_cd_from_histo(l3D0Rate_motherBin_2,l3D0Rate_motherBin_2_cd);
  make_cd_from_histo(l3D0Rate_motherBin_3,l3D0Rate_motherBin_3_cd);
  make_cd_from_histo(l3D0Rate_motherBin_4,l3D0Rate_motherBin_4_cd);
  make_cd_from_histo(l3D0Rate_motherBin_5,l3D0Rate_motherBin_5_cd);
  make_cd_from_histo(l3D0Rate_motherBin_6,l3D0Rate_motherBin_6_cd);
  make_cd_from_histo(l3D0Rate_motherBin_10,l3D0Rate_motherBin_10_cd);
  make_cd_from_histo(l3D0Rate_motherBin_11,l3D0Rate_motherBin_11_cd);
  make_cd_from_histo(l3D0Rate_motherBin_12,l3D0Rate_motherBin_12_cd);
  make_cd_from_histo(l3D0Rate_motherBin_13,l3D0Rate_motherBin_13_cd);
  make_cd_from_histo(l3D0Rate_motherBin_14,l3D0Rate_motherBin_14_cd);
  make_cd_from_histo(l3D0Rate_motherBin_15,l3D0Rate_motherBin_15_cd);

  if(!isData){
    l3D0Rate->Add(l3D0Rate_motherBin_1);
    l3D0Rate->Add(l3D0Rate_motherBin_2);
    //l3D0Rate->Add(l3D0Rate_motherBin_3);
    //l3D0Rate->Add(l3D0Rate_motherBin_4);
    //l3D0Rate->Add(l3D0Rate_motherBin_5);
    //l3D0Rate->Add(l3D0Rate_motherBin_6);
    l3D0Rate->Add(l3D0Rate_motherBin_10);
    l3D0Rate->Add(l3D0Rate_motherBin_11);
    l3D0Rate->Add(l3D0Rate_motherBin_12);
    l3D0Rate->Add(l3D0Rate_motherBin_13);
    l3D0Rate->Add(l3D0Rate_motherBin_14);
  } else {
    l3D0Rate->Add(l3D0Rate_motherBin_15);
  }

  if(!isData){
    l3D0Rate_cd->Add(l3D0Rate_motherBin_1_cd);
    l3D0Rate_cd->Add(l3D0Rate_motherBin_2_cd);
    //l3D0Rate_cd->Add(l3D0Rate_motherBin_3_cd);
    //l3D0Rate_cd->Add(l3D0Rate_motherBin_4_cd);
    //l3D0Rate_cd->Add(l3D0Rate_motherBin_5_cd);
    //l3D0Rate_cd->Add(l3D0Rate_motherBin_6_cd);
    l3D0Rate_cd->Add(l3D0Rate_motherBin_10_cd);
    l3D0Rate_cd->Add(l3D0Rate_motherBin_11_cd);
    l3D0Rate_cd->Add(l3D0Rate_motherBin_12_cd);
    l3D0Rate_cd->Add(l3D0Rate_motherBin_13_cd);
    l3D0Rate_cd->Add(l3D0Rate_motherBin_14_cd);
  } else {
    l3D0Rate_cd->Add(l3D0Rate_motherBin_15_cd);
  }
  //aaa
  
  fill_overflow(l2D0Rate_motherBin_1,l2OverFlow[0]);
  fill_overflow(l2D0Rate_motherBin_2,l2OverFlow[1]);
  fill_overflow(l2D0Rate_motherBin_3,l2OverFlow[2]);
  fill_overflow(l2D0Rate_motherBin_4,l2OverFlow[3]);
  fill_overflow(l2D0Rate_motherBin_5,l2OverFlow[4]);
  fill_overflow(l2D0Rate_motherBin_6,l2OverFlow[5]);
  fill_overflow(l2D0Rate_motherBin_10,l2OverFlow[9]);
  fill_overflow(l2D0Rate_motherBin_11,l2OverFlow[10]);
  fill_overflow(l2D0Rate_motherBin_12,l2OverFlow[11]);
  fill_overflow(l2D0Rate_motherBin_13,l2OverFlow[12]);
  fill_overflow(l2D0Rate_motherBin_14,l2OverFlow[13]);
  fill_overflow(l2D0Rate_motherBin_15,l2OverFlow[14]);

  make_cd_from_histo(l2D0Rate_motherBin_1,l2D0Rate_motherBin_1_cd);
  make_cd_from_histo(l2D0Rate_motherBin_2,l2D0Rate_motherBin_2_cd);
  make_cd_from_histo(l2D0Rate_motherBin_3,l2D0Rate_motherBin_3_cd);
  make_cd_from_histo(l2D0Rate_motherBin_4,l2D0Rate_motherBin_4_cd);
  make_cd_from_histo(l2D0Rate_motherBin_5,l2D0Rate_motherBin_5_cd);
  make_cd_from_histo(l2D0Rate_motherBin_6,l2D0Rate_motherBin_6_cd);
  make_cd_from_histo(l2D0Rate_motherBin_10,l2D0Rate_motherBin_10_cd);
  make_cd_from_histo(l2D0Rate_motherBin_11,l2D0Rate_motherBin_11_cd);
  make_cd_from_histo(l2D0Rate_motherBin_12,l2D0Rate_motherBin_12_cd);
  make_cd_from_histo(l2D0Rate_motherBin_13,l2D0Rate_motherBin_13_cd);
  make_cd_from_histo(l2D0Rate_motherBin_14,l2D0Rate_motherBin_14_cd);
  make_cd_from_histo(l2D0Rate_motherBin_15,l2D0Rate_motherBin_15_cd);

  if(!isData){
  l2D0Rate->Add(l2D0Rate_motherBin_1);
  l2D0Rate->Add(l2D0Rate_motherBin_2);
  //l2D0Rate->Add(l2D0Rate_motherBin_3);
  //l2D0Rate->Add(l2D0Rate_motherBin_4);
  //l2D0Rate->Add(l2D0Rate_motherBin_5);
  //l2D0Rate->Add(l2D0Rate_motherBin_6);
  l2D0Rate->Add(l2D0Rate_motherBin_10);
  l2D0Rate->Add(l2D0Rate_motherBin_11);
  l2D0Rate->Add(l2D0Rate_motherBin_12);
  l2D0Rate->Add(l2D0Rate_motherBin_13);
  l2D0Rate->Add(l2D0Rate_motherBin_14);
  } else {
    l2D0Rate->Add(l2D0Rate_motherBin_15);
  }

  if(!isData){
    l2D0Rate_cd->Add(l2D0Rate_motherBin_1_cd);
    l2D0Rate_cd->Add(l2D0Rate_motherBin_2_cd);
    //l2D0Rate_cd->Add(l2D0Rate_motherBin_3_cd);
    //l2D0Rate_cd->Add(l2D0Rate_motherBin_4_cd);
    //l2D0Rate_cd->Add(l2D0Rate_motherBin_5_cd);
    //l2D0Rate_cd->Add(l2D0Rate_motherBin_6_cd);
    l2D0Rate_cd->Add(l2D0Rate_motherBin_10_cd);
    l2D0Rate_cd->Add(l2D0Rate_motherBin_11_cd);
    l2D0Rate_cd->Add(l2D0Rate_motherBin_12_cd);
    l2D0Rate_cd->Add(l2D0Rate_motherBin_13_cd);
    l2D0Rate_cd->Add(l2D0Rate_motherBin_14_cd);
  } else {
    l2D0Rate_cd->Add(l2D0Rate_motherBin_15_cd);
  }
  //aaa

  fill_overflow(tkTrackD0Rate_motherBin_1,tkTrackOverFlow[0]);
  fill_overflow(tkTrackD0Rate_motherBin_2,tkTrackOverFlow[1]);
  fill_overflow(tkTrackD0Rate_motherBin_3,tkTrackOverFlow[2]);
  fill_overflow(tkTrackD0Rate_motherBin_4,tkTrackOverFlow[3]);
  fill_overflow(tkTrackD0Rate_motherBin_5,tkTrackOverFlow[4]);
  fill_overflow(tkTrackD0Rate_motherBin_6,tkTrackOverFlow[5]);
  fill_overflow(tkTrackD0Rate_motherBin_10,tkTrackOverFlow[9]);
  fill_overflow(tkTrackD0Rate_motherBin_11,tkTrackOverFlow[10]);
  fill_overflow(tkTrackD0Rate_motherBin_12,tkTrackOverFlow[11]);
  fill_overflow(tkTrackD0Rate_motherBin_13,tkTrackOverFlow[12]);
  fill_overflow(tkTrackD0Rate_motherBin_14,tkTrackOverFlow[13]);
  fill_overflow(tkTrackD0Rate_motherBin_15,tkTrackOverFlow[14]);

  make_cd_from_histo(tkTrackD0Rate_motherBin_1,tkTrackD0Rate_motherBin_1_cd);
  make_cd_from_histo(tkTrackD0Rate_motherBin_2,tkTrackD0Rate_motherBin_2_cd);
  make_cd_from_histo(tkTrackD0Rate_motherBin_3,tkTrackD0Rate_motherBin_3_cd);
  make_cd_from_histo(tkTrackD0Rate_motherBin_4,tkTrackD0Rate_motherBin_4_cd);
  make_cd_from_histo(tkTrackD0Rate_motherBin_5,tkTrackD0Rate_motherBin_5_cd);
  make_cd_from_histo(tkTrackD0Rate_motherBin_6,tkTrackD0Rate_motherBin_6_cd);
  make_cd_from_histo(tkTrackD0Rate_motherBin_10,tkTrackD0Rate_motherBin_10_cd);
  make_cd_from_histo(tkTrackD0Rate_motherBin_11,tkTrackD0Rate_motherBin_11_cd);
  make_cd_from_histo(tkTrackD0Rate_motherBin_12,tkTrackD0Rate_motherBin_12_cd);
  make_cd_from_histo(tkTrackD0Rate_motherBin_13,tkTrackD0Rate_motherBin_13_cd);
  make_cd_from_histo(tkTrackD0Rate_motherBin_14,tkTrackD0Rate_motherBin_14_cd);
  make_cd_from_histo(tkTrackD0Rate_motherBin_15,tkTrackD0Rate_motherBin_15_cd);

  if(!isData){
    tkTrackD0Rate->Add(tkTrackD0Rate_motherBin_1);
    tkTrackD0Rate->Add(tkTrackD0Rate_motherBin_2);
    //tkTrackD0Rate->Add(tkTrackD0Rate_motherBin_3);
    //tkTrackD0Rate->Add(tkTrackD0Rate_motherBin_4);
    //tkTrackD0Rate->Add(tkTrackD0Rate_motherBin_5);
    //tkTrackD0Rate->Add(tkTrackD0Rate_motherBin_6);
    tkTrackD0Rate->Add(tkTrackD0Rate_motherBin_10);
    tkTrackD0Rate->Add(tkTrackD0Rate_motherBin_11);
    tkTrackD0Rate->Add(tkTrackD0Rate_motherBin_12);
    tkTrackD0Rate->Add(tkTrackD0Rate_motherBin_13);
    tkTrackD0Rate->Add(tkTrackD0Rate_motherBin_14);
  } else {
    tkTrackD0Rate->Add(tkTrackD0Rate_motherBin_15);
  }
  
  if(!isData){
    tkTrackD0Rate_cd->Add(tkTrackD0Rate_motherBin_1_cd);
    tkTrackD0Rate_cd->Add(tkTrackD0Rate_motherBin_2_cd);
    //tkTrackD0Rate_cd->Add(tkTrackD0Rate_motherBin_3_cd);
    //tkTrackD0Rate_cd->Add(tkTrackD0Rate_motherBin_4_cd);
    //tkTrackD0Rate_cd->Add(tkTrackD0Rate_motherBin_5_cd);
    //tkTrackD0Rate_cd->Add(tkTrackD0Rate_motherBin_6_cd);
    tkTrackD0Rate_cd->Add(tkTrackD0Rate_motherBin_10_cd);
    tkTrackD0Rate_cd->Add(tkTrackD0Rate_motherBin_11_cd);
    tkTrackD0Rate_cd->Add(tkTrackD0Rate_motherBin_12_cd);
    tkTrackD0Rate_cd->Add(tkTrackD0Rate_motherBin_13_cd);
    tkTrackD0Rate_cd->Add(tkTrackD0Rate_motherBin_14_cd);
  } else {
    tkTrackD0Rate_cd->Add(tkTrackD0Rate_motherBin_15_cd);
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
  scaleToRate("*D0Rate*", rateFactor);


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
