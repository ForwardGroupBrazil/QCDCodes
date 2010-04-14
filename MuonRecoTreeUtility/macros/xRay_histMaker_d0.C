#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
//#include "Math/LorentzVector.h"
#include <vector>
#include <map>
#include "RecoTree.h"
#include <stdio.h>

typedef std::pair< std::string, std::vector<int> > bin;
std::vector<bin> Bins;
//in which bin a pdgID goes
std::map<int, int> mymap;

std::map<std::string,TH1*> muonPtHistoMap;
std::map<std::string,TH1*> l2PtHistoMap;
std::map<std::string,TH1*> tkPtHistoMap;

//int ScanTree ( TTree* tree) {
int ScanTreeD0 ( TTree* tree, char *fileName, bool isData=false,double weight = 1.0) {
  
  initArrays();
  
  // This reads in the tree.  As you might imagine.
  Init(tree);
  TFile *histFile = new TFile(fileName,"UPDATE");
  TDirectory *histDir = histFile->mkdir("D0Rate");
  
  histDir->cd();
  
  // Time for some histograms and stacks
  
  double ptLowerLimit = -5.;
  double ptUpperLimit = 5.;
  int nBins = 50;

  //Book the L3 histos  
  for (int i = 0; i != Bins.size(); ++i) {    
    char histoName[20];
    sprintf(histoName,"l3D0Rate_%d",i);
    TString histoString(histoName);
    // cout << "histoName : " << histoName << endl;
    // cout << "histoString : " << histoString.Data() << endl;
    TString histoTitle(GetBinName(i));
    TH1 * tmpTh1 = new TH1F(histoString,histoTitle,nBins,ptLowerLimit,ptUpperLimit);
    muonPtHistoMap[histoName] = tmpTh1;
  }
  
  char histoName[20];
  sprintf(histoName,"l3D0Rate_%d",Bins.size());
  TString histoString(histoName);
  // cout << "histoName : " << histoName << endl;
  // cout << "histoString : " << histoString.Data() << endl;
  TString histoTitle("Overflow");
  TH1 * tmpTh1 = new TH1F(histoString,histoTitle,nBins,ptLowerLimit,ptUpperLimit);
  muonPtHistoMap[histoName] = tmpTh1;
  
  for (int i = 0; i != Bins.size()+1; ++i) {
    char histoName[20];
    sprintf(histoName,"l3D0Rate_%d",i);
    cout << "Histo " << i << " " << muonPtHistoMap[histoName]->GetName() << " " << muonPtHistoMap[histoName]->GetTitle() << endl;
  }

  //Book the L2 histos  
  for (int i = 0; i != Bins.size(); ++i) {    
    char histoName[20];
    sprintf(histoName,"l2D0Rate_%d",i);
    TString histoString(histoName);
    // cout << "histoName : " << histoName << endl;
    // cout << "histoString : " << histoString.Data() << endl;
    TString histoTitle(GetBinName(i));
    TH1 * tmpTh1 = new TH1F(histoString,histoTitle,nBins,ptLowerLimit,ptUpperLimit);
    l2PtHistoMap[histoName] = tmpTh1;
  }
  
  char histoName2[20];
  sprintf(histoName2,"l2D0Rate_%d",Bins.size());
  TString histoString2(histoName2);
  // cout << "histoName : " << histoName2 << endl;
  // cout << "histoString : " << histoString2.Data() << endl;
  TString histoTitle2("Overflow");
  TH1 * tmpTh2 = new TH1F(histoString2,histoTitle2,nBins,ptLowerLimit,ptUpperLimit);
  l2PtHistoMap[histoName2] = tmpTh2;
  
  for (int i = 0; i != Bins.size()+1; ++i) {
    char histoName[20];
    sprintf(histoName,"l2D0Rate_%d",i);
    cout << "Histo " << i << " " << l2PtHistoMap[histoName]->GetName() << " " << l2PtHistoMap[histoName]->GetTitle() << endl;
  }

  //Book the Tk histos  
  for (int i = 0; i != Bins.size(); ++i) {    
    char histoName[20];
    sprintf(histoName,"tkTrackD0Rate_%d",i);
    TString histoString(histoName);
    // cout << "histoName : " << histoName << endl;
    // cout << "histoString : " << histoString.Data() << endl;
    TString histoTitle(GetBinName(i));
    TH1 * tmpTh1 = new TH1F(histoString,histoTitle,nBins,ptLowerLimit,ptUpperLimit);
    tkPtHistoMap[histoName] = tmpTh1;
  }
  
  char histoName3[20];
  sprintf(histoName3,"tkTrackD0Rate_%d",Bins.size());
  TString histoString3(histoName3);
  // cout << "histoName : " << histoName3 << endl;
  // cout << "histoString : " << histoString3.Data() << endl;
  TString histoTitle3("Overflow");
  TH1 * tmpTh3 = new TH1F(histoString3,histoTitle3,nBins,ptLowerLimit,ptUpperLimit);
  tkPtHistoMap[histoName3] = tmpTh3;
  
  for (int i = 0; i != Bins.size()+1; ++i) {
    char histoName[20];
    sprintf(histoName,"tkTrackD0Rate_%d",i);
    cout << "Histo " << i << " " << tkPtHistoMap[histoName]->GetName() << " " << tkPtHistoMap[histoName]->GetTitle() << endl;
  }
    
  THStack *l3D0Rate = new THStack("l3D0Rate","L3 rate as f(p_{T,L3})");
  //THStack *l3D0Rate_cd = new THStack("l3D0Rate_cd","L3 rate as f(p_{T,L3})");
    
  THStack *l2D0Rate = new THStack("l2D0Rate","L2 rate as f(p_{T,L2})");
  //THStack *l2D0Rate_cd = new THStack("l2D0Rate_cd","L2 rate as f(p_{T,L2})");
  
  THStack *tkTrackD0Rate = new THStack("tkTrackD0Rate","tkTrack rate as f(p_{T,tkTrack})");
  //THStack *tkTrackD0Rate_cd = new THStack("tkTrackD0Rate_cd","tkTrack rate as f(p_{T,tkTrack})");
  
  // There are many comments to be made about ROOT.  Most of them Rated-R or Rated-X.
  // Here, I'll limit myself to saying this is necessary for the TStacks to go into the TDirectory.
  gDirectory->Append(l2D0Rate);
  gDirectory->Append(l3D0Rate);
  gDirectory->Append(tkTrackD0Rate);
  
  //gDirectory->Append(l2D0Rate_cd);
  //gDirectory->Append(l3D0Rate_cd);
  //gDirectory->Append(tkTrackD0Rate_cd);
  
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

  FILE * pFile;  
  pFile = fopen ("myL3ParentIdBin.txt","w");
  
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
      
      //if ((*muTrackerMuonArbitrated).at(iMu)) passTK = true;
      
      if (   (*muAllGlobalMuons).at(iMu)
      	     ) {passL3 = true; passTK = true;}

      int myIsAssociated = (*l3IsAssociated).at(iMu);
      int myL3AssociationPdgId = (*l3AssociationPdgId).at(iMu);
      int myL3ParentID = (*l3ParentID).at(iMu);
      int myL3MotherBinNumber = (*l3MotherBinNumber).at(iMu);
      int myL3TPIdBin = GetBinNum((*l3AssociationPdgId).at(iMu));
      int myL3ParentIdBin = GetBinNum((*l3ParentID).at(iMu));

      /*
      if (passL3) {
	cout << "iEntry " << iEntry 
	     << " iMu " << iMu  << " of " << nMu 
	     << " isAssoc " << myIsAssociated 
	     << " pdgId " << myL3AssociationPdgId
	     << " parentID " << myL3ParentID 
	     << " mother bin " << myL3MotherBinNumber
	     << endl;
      }
      */

      if (passL3) {
	fprintf(pFile,"i %d iMu %d of %d isAssoc %d pdgId %d parentID %d mother bin %d lBin %d %d\n", iEntry, iMu, nMu, myIsAssociated, myL3AssociationPdgId, myL3ParentID, myL3MotherBinNumber,myL3TPIdBin,myL3ParentIdBin);
      }
      
      /*
	if(myL3TPIdBin == 10) {
	fprintf(pFile,"xEntry %d iMu %d of %d isAssoc %d pdgId %d parentID %d mother bin %d \n", iEntry, iMu, nMu, myIsAssociated, myL3AssociationPdgId, myL3ParentID, myL3MotherBinNumber);
	}
	
	if(myL3ParentIdBin == 10) {
	fprintf(pFile,"yEntry %d iMu %d of %d isAssoc %d pdgId %d parentID %d mother bin %d \n", iEntry, iMu, nMu, myIsAssociated, myL3AssociationPdgId, myL3ParentID, myL3MotherBinNumber);
	}
      */
      
      if (passL3) {
	//double pt_L3 = findMaxPt(l3Pt,l3Eta,l3D0);
	//int l3_index = findIndexOfMaxPt(l3Pt,l3Eta,l3D0);
	//if(!(*muAllArbitrated).at(iMu)) break;
	int l3_index = iMu;
	double pt_L3 = (*l3D0).at(iMu);
		
	if(!isData) {
	  //start adam test
	  if( ( myIsAssociated == 1 ) && 
	      (
	       myL3MotherBinNumber==777 ||
	       myL3MotherBinNumber==-777 ||
	       //myL3MotherBinNumber==-888 ||
	       myL3MotherBinNumber==-999	
	       )
	      ) {	    
	    fprintf(pFile,"j %d iMu %d of %d isAssoc %d pdgId %d parentID %d mother bin %d lBin %d %d\n", iEntry, iMu, nMu, myIsAssociated, myL3AssociationPdgId, myL3ParentID, myL3MotherBinNumber,myL3TPIdBin,myL3ParentIdBin);	    
	  }
	  
	  // end adam test
	  
	  //int myL3TPId = myL3TPIdBin;
	  //int myL3ParentId = myL3ParentIdBin;
	  //int myMotherBinNumber = myL3MotherBinNumber;

	  //adam test 2
	  //if( fabs(myL3AssociationPdgId) != 13) {
	  char histoName[20];
	  if(myL3ParentIdBin==Bins.size()) myL3ParentIdBin = myL3ParentIdBin-2;
	  sprintf(histoName,"l3D0Rate_%d",myL3ParentIdBin);	  
	  muonPtHistoMap[histoName]->Fill(pt_L3,1.0*weight);
	  //}
	  //adam test 2

	  /*
	  if (myIsAssociated == 1 
	      && (*l3AssociationVar).at(l3_index) > -0.1 
	      && fabs(myL3AssociationPdgId) == 13 ) {
	    
	    char histoName[20];
	    sprintf(histoName,"l3D0Rate_%d",myL3ParentIdBin);
	    
	    muonPtHistoMap[histoName]->Fill(pt_L3,weight);
	    
	    if(
	       myL3MotherBinNumber==777 ||
	       myL3MotherBinNumber==-777 ||
	       myL3MotherBinNumber==-888 ||
	       myL3MotherBinNumber==-999	       
	       ) {
	      fprintf(pFile,"k %d iMu %d of %d isAssoc %d pdgId %d parentID %d mother bin %d lBin %d %d\n", iEntry, iMu, nMu, myIsAssociated, myL3AssociationPdgId, myL3ParentID, myL3MotherBinNumber,myL3TPIdBin,myL3ParentIdBin);
	    } 
	  } else if (myIsAssociated == 1 && (*l3AssociationVar).at(l3_index) > -0.1 && fabs(myL3AssociationPdgId)!=13 && myL3MotherBinNumber == -888) {	    
	    //muonPtHistoMap["l3D0Rate_9"]->Fill(pt_L3,weight);
	    char histoName[20];
	    sprintf(histoName,"l3D0Rate_%d",myL3ParentIdBin);
	    
	    muonPtHistoMap[histoName]->Fill(pt_L3,weight);
	  } else if (myIsAssociated == 1 && (*l3AssociationVar).at(l3_index) > -0.1 && fabs(myL3AssociationPdgId)!=13 && myL3MotherBinNumber == -777) {
	    muonPtHistoMap["l3D0Rate_8"]->Fill(pt_L3,weight);
	  } else {
	    muonPtHistoMap["l3D0Rate_9"]->Fill(pt_L3,weight);
	  } */
	} else {
	  //cout << "line 273 Data" << endl;
	  char histoName[20];
	  //if(myL3ParentIdBin==Bins.size()) myL3ParentIdBin--;
	  sprintf(histoName,"l3D0Rate_%d",Bins.size()-1);	  
	  muonPtHistoMap[histoName]->Fill(pt_L3,weight);	  
	  //muonPtHistoMap["l3D0Rate_9"]->Fill(pt_L3,weight);
	}
      }
      
      if (passL2) {
	//double pt_L2 = findMaxPt(l2Pt,l2Eta,l2D0);
	//int l2_index = findIndexOfMaxPt(l2Pt,l2Eta,l2D0);
	int l2_index = iMu;
	double pt_L2 = (*l2D0).at(iMu);
	if(!isData) {

      }
      }

      if (passTK) {
	//double pt_tkTrack = findMaxPt(tkTrackPt,tkTrackEta,tkTrackD0);
	//int tkTrack_index = findIndexOfMaxPt(tkTrackPt,tkTrackEta,tkTrackD0);
	int tkTrack_index = iMu;
	double pt_tkTrack = (*tkTrackD0).at(iMu);

	int myIsAssociated = (*tkTrackIsAssociated).at(iMu);
	int myL3AssociationPdgId = myIsAssociated ? (*tkTrackAssociationPdgId).at(iMu) : -777;
	int myL3ParentID = (*tkTrackParentID).at(iMu);
	int myL3MotherBinNumber = (*tkTrackMotherBinNumber).at(iMu);
	int myL3TPIdBin = myIsAssociated ? GetBinNum((*tkTrackAssociationPdgId).at(iMu)) : -777;
	int myL3ParentIdBin = GetBinNum((*tkTrackParentID).at(iMu));


	if(!isData) {
	  char histoName[20];
	  if(myL3ParentIdBin==Bins.size()) myL3ParentIdBin = myL3ParentIdBin-2;
	  sprintf(histoName,"tkTrackD0Rate_%d",myL3ParentIdBin);	  
	  tkPtHistoMap[histoName]->Fill(pt_tkTrack,1.0*weight);	
	} else {
	  //cout << "line 273 Data" << endl;
	  char histoName[20];
	  //if(myL3ParentIdBin==Bins.size()) myL3ParentIdBin--;
	  sprintf(histoName,"tkTrackD0Rate_%d",Bins.size()-1);	  
	  tkPtHistoMap[histoName]->Fill(pt_tkTrack,weight);	  
	  //muonPtHistoMap["l3D0Rate_9"]->Fill(pt_L3,weight);
	}
      }
    }//cout << "line 399 loop over muons" << endl;
  }cout << "line 400 loop over entries" << endl;
  
     
  if(!isData){
    //for (int i = 0; i != Bins.size()+1; ++i) {
    for (int i = 1; i != Bins.size()-1; ++i) {
      char histoName[20];
      sprintf(histoName,"l3D0Rate_%d",i);
      l3D0Rate->Add(muonPtHistoMap[histoName]);
    }
  } else {
    char histoName[20];
    sprintf(histoName,"l3D0Rate_%d",Bins.size()-1);
    l3D0Rate->Add(muonPtHistoMap[histoName]);
  }

  if(!isData){
    //for (int i = 0; i != Bins.size()+1; ++i) {
    for (int i = 1; i != Bins.size()-1; ++i) {
      char histoName[20];
      sprintf(histoName,"l2D0Rate_%d",i);
      l2D0Rate->Add(l2PtHistoMap[histoName]);
    }
  } else {
    char histoName[20];
    sprintf(histoName,"l2D0Rate_%d",Bins.size()-1);
    l2D0Rate->Add(l2PtHistoMap[histoName]);
  }

  if(!isData){
    //for (int i = 0; i != Bins.size()+1; ++i) {
    for (int i = 1; i != Bins.size()-1; ++i) {
      char histoName[20];
      sprintf(histoName,"tkTrackD0Rate_%d",i);
      tkTrackD0Rate->Add(tkPtHistoMap[histoName]);
    }
  } else {
    char histoName[20];
    sprintf(histoName,"tkTrackD0Rate_%d",Bins.size()-1);
    tkTrackD0Rate->Add(tkPtHistoMap[histoName]);
  }

  double crossSection_mb = 51.6; // ppMuX
  double filterEff = 0.00289; // ppMuX

  double eventNumberUnit= 1000000;

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
  //scaleToRate("*D0Rate*", rateFactor);


  histDir->Write("",TObject::kOverwrite);

  histFile->Close();

  fclose (pFile);
  
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
    //aaa cout << "Testing name: " << name << " against " << collectionName << endl;
    
    if (TString(collectionName).MaybeRegexp()) {
      //aaa cout << "we have a possible match" << endl;
      //aaa cout << "Trying to match to " << TString(obj->GetName()) << endl;
      if (TString(obj->GetName()).Index(reg) < 0 ) {
	cout << "failure here.  Argument returns " << TString(obj->GetName()).Index(reg) << endl;
	continue;
      }
    }
    else if (! name.BeginsWith(collectionName)) continue;
    
    //aaa cout << "We're trying to scale" << name << endl;
    ((TH1*)obj)->Scale(rateFactor);
    
  }
  
}

void initArrays()
{
  //configured by hand as Finn did it
  Bins.resize(10);

  Bins[0].first = "ID=0";
  Bins[0].second.push_back(0);

  Bins[1].first = "#pi+/-";
  Bins[1].second.push_back(211);
  Bins[1].second.push_back(-211);

  Bins[2].first = "K";
  Bins[2].second.push_back(321);
  Bins[2].second.push_back(-321);
  Bins[2].second.push_back(130);
  Bins[2].second.push_back(-130);

  Bins[3].first = "D";
  Bins[3].second.push_back(411);
  Bins[3].second.push_back(-411);
  Bins[3].second.push_back(421);
  Bins[3].second.push_back(-421);
  Bins[3].second.push_back(431);
  Bins[3].second.push_back(-431);

  Bins[4].first = "B";
  Bins[4].second.push_back(521);
  Bins[4].second.push_back(-521);
  Bins[4].second.push_back(511);
  Bins[4].second.push_back(-511);
  Bins[4].second.push_back(531);
  Bins[4].second.push_back(-531);

  Bins[5].first = "#Lambda_{b}";
  Bins[5].second.push_back(5122);
  Bins[5].second.push_back(-5122);

  Bins[6].first = "J/#Psi";
  Bins[6].second.push_back(443);
  Bins[6].second.push_back(-443);

  Bins[7].first = "#Upsilon(nS)";
  Bins[7].second.push_back(553);
  Bins[7].second.push_back(-553);
  Bins[7].second.push_back(200553);
  Bins[7].second.push_back(-200553);
  Bins[7].second.push_back(100553);
  Bins[7].second.push_back(-100553);
  Bins[7].second.push_back(300553);
  Bins[7].second.push_back(-300553);
  
  //Bins[8].first = "W^{+/-}";
  //Bins[8].second.push_back(24);
  //Bins[8].second.push_back(-24);
  
  //Bins[9].first = "Z^{0}";
  //Bins[9].second.push_back(23);
  
  //Bins[10].first = "#tau";
  //Bins[10].second.push_back(15);
  //Bins[10].second.push_back(-15);
  
  Bins[8].first = "Fakes";
  Bins[8].second.push_back(-777);

  Bins[9].first = "Data";
  Bins[9].second.push_back(999);
  
  //Bins[9].first = "-888";
  //Bins[9].second.push_back(-888);
  
  //Bins[10].first = "-999";
  //Bins[10].second.push_back(-999);
  
  //Bins[11].first = "+777";
  //Bins[11].second.push_back(777);
  
  //and filling the map
  for (int r=0;r!=Bins.size();++r){
    for (int i=0;i!=Bins[r].second.size();++i)
      {mymap[Bins[r].second[i]] = r;}
  }
  //mymap[-999] = Bins.size();
}

int GetBinNum(int pdgID)
{
  std::map<int, int>::iterator i=mymap.find(pdgID);
  //if registerd pdgid, return the bin
  if ( i!= mymap.end()){return i->second;}
  //else return overflow bin
  else {
    //cout << "IDconverttoBinNum " << pdgID << " is not a registered pdgID." << endl;
    return Bins.size();}
}

std::string GetBinName(int BinNumvalue)
{
  if (BinNumvalue>Bins.size()) return "bin number is not valid";
  else if (BinNumvalue==Bins.size()) return other();
  else { return Bins[BinNumvalue].first; }
}
