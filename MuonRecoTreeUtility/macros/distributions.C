#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
//#include "Math/LorentzVector.h"
#include <vector.h>
#include "RecoTree.h"
#include <map.h>

#include <sstream>

template <class T>
inline std::string to_string (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}

int distributions ( TTree* tree, char *fileName="histo.root", double effCS=1.0) {

  TFile *histFileAlgo = new TFile(fileName,"RECREATE");
  TDirectory *histDirAlgo = histFileAlgo->mkdir("test_hist");

  Init(tree);

  TH1F* glb_pt[8];
  TH1F* glb_deltaPt[8];
  TH1F* glb_calComp[8];
  TH1F* glb_segComp[8];
  TH1F* glb_nChamber[8];
  TH1F* glb_nChamberMatch[8];
  TH1F* glb_NCSC_1[8];
  TH1F* glb_NCSC_2[8];
  TH1F* glb_NCSC_3[8];
  TH1F* glb_NCSC_4[8];
  TH1F* glb_NDT_1[8];
  TH1F* glb_NDT_2[8];
  TH1F* glb_NDT_3[8];
  TH1F* glb_NDT_4[8];
  TH1F* glb_NCSC_total[8];
  TH1F* glb_NDT_total[8];
  TH1F* glb_NSeg_total[8];
  TH1F* glb_muTrkKink[8];
  TH1F* glb_muGlbKink[8];
  TH1F* glb_muTrkRelChi2[8];
  TH1F* glb_muStaRelChi2[8];
  TH1F* glb_muDeltaTrkRelChi2[8];
  TH1F* glb_muDeltaStaRelChi2[8];
  TH1F* glb_nChi2[8];

  TH1F* mu_TMLastStationLoose[8];
  TH1F* mu_TMLastStationTight[8];
  TH1F* mu_TM2DCompatibilityLoose[8];
  TH1F* mu_TM2DCompatibilityTight[8];
  TH1F* mu_GlobalMuonPromptTight[8];

  TDirectory* dirs[8];
  TH1F* tmp;
  histDirAlgo->cd();
  for (unsigned int ww=0;ww<8;ww++){
    histDirAlgo->cd();
    std::string ss = to_string(ww);

    TString name("dir");
    name+=ss;
    TString name2("l3Pt");
    //name2+=ss;
    TDirectory *tmpDir = histDirAlgo->mkdir(name.Data());
    dirs[ww] = tmpDir;
    tmpDir->cd();
    tmp = new TH1F(name2,"histogram of L3 p_{T}",250,0,500);
    tmp->SetDirectory(tmpDir);
    glb_pt[ww] = tmp;
    glb_deltaPt[ww] = new TH1F("deltaPt","histogram of L3 p_{T}",100,-10,10);
    mu_TMLastStationLoose[ww] = new TH1F("mu_TMLastStationLoose","mu_TMLastStationLoose",3,-1.5,1.5);
    mu_TMLastStationTight[ww] = new TH1F("mu_TMLastStationTight","mu_TMLastStationTight",3,-1.5,1.5);;
    mu_TM2DCompatibilityLoose[ww] = new TH1F("mu_TM2DCompatibilityLoose","mu_TM2DCompatibilityLoose",3,-1.5,1.5);;
    mu_TM2DCompatibilityTight[ww] = new TH1F("mu_TM2DCompatibilityTight","mu_TM2DCompatibilityTight",3,-1.5,1.5);;
    mu_GlobalMuonPromptTight[ww] = new TH1F("mu_GlobalMuonPromptTight","mu_GlobalMuonPromptTight",3,-1.5,1.5);;
    glb_calComp[ww] = new TH1F("glb_calComp","Calo Compatibility",20,0.,1.);
    glb_segComp[ww] = new TH1F("glb_segComp","Segment Compatibility",20,0.,1.);
    glb_nChamber[ww] = new TH1F("glb_nChamber","Number of Chambers",11,-0.5,10.5);
    glb_nChamberMatch[ww] = new TH1F("glb_nChamberMatch","Number of Matched Chambers",11,-0.5,10.5);
    //glb_stationMask[ww] = new TH1F("glb_stationMask","Station Mask",11,-0.5,10.5);
    glb_NCSC_1[ww] = new TH1F("glb_nCSC_1","Number of segments in station 1",11,-0.5,10.5);
    glb_NCSC_2[ww] = new TH1F("glb_nCSC_2","Number of segments in station 2",11,-0.5,10.5);
    glb_NCSC_3[ww] = new TH1F("glb_nCSC_3","Number of segments in station 3",11,-0.5,10.5);
    glb_NCSC_4[ww] = new TH1F("glb_nCSC_4","Number of segments in station 4",11,-0.5,10.5);
    glb_NDT_1[ww] = new TH1F("glb_nDT_1","Number of segments in DT station 1",11,-0.5,10.5);
    glb_NDT_2[ww] = new TH1F("glb_nDT_2","Number of segments in DT station 2",11,-0.5,10.5);
    glb_NDT_3[ww] = new TH1F("glb_nDT_3","Number of segments in DT station 3",11,-0.5,10.5);
    glb_NDT_4[ww] = new TH1F("glb_nDT_4","Number of segments in DT station 4",11,-0.5,10.5);
    glb_NCSC_total[ww] = new TH1F("glb_nCSC_total","Number of segments in all CSC stations",11,-0.5,10.5);
    glb_NDT_total[ww] = new TH1F("glb_nDT_total","Number of segments in all DT stations",11,-0.5,10.5);
    glb_NSeg_total[ww] = new TH1F("glb_nSeg_total","Number of segments in all stations",16,-0.5,15.5);
    glb_muTrkKink[ww] = new TH1F("glb_trkKink","Tracker Kink",500,0,500);
    glb_muGlbKink[ww] = new TH1F("glb_glbKink","Global Kink",500,0,500);
    glb_muTrkRelChi2[ww] = new TH1F("glb_trkRelChi2","Tracker Relative #Chi^{2}",50,0,50);
    glb_muStaRelChi2[ww] = new TH1F("glb_staRelChi2","Standalone Relative #Chi^{2}",50,0,50);
    glb_muDeltaTrkRelChi2[ww] = new TH1F("glb_trkDeltaRelChi2","#Delta Tracker Relative #Chi^{2}",10,0,5);
    glb_muDeltaStaRelChi2[ww] = new TH1F("glb_staDeltaRelChi2","#Delta Standalone Relative #Chi^{2}",10,0,5);
    glb_nChi2[ww] = new TH1F("glb_nChi2","Global #Chi^{2}/DoF",50,0,50);
  }
  histDirAlgo->cd();

  TH1F *l3TestPtAll = new TH1F("l3TestPtAll","histogram of L3 p_{T}",100,0,100);
  TH1F *l3deltaPt = new TH1F("l3deltaPt","histogram of L3 #Delta p_{T}",100,-10.,10.);

  //TH1F *recHitsTest = new TH1F("recHitsTest","test to fill x position of recHits",1600,-800,800); 
  //TH3F *poorMansDisplay = new TH3F("poorMansDisplay","test for hits", 200, -800.0, 800.0, 200, -800.0, 800.0, 250, -1045.0, 1045.0); 

  int nEntries = tree->GetEntries();
  double weight = nEntries > 0 ? effCS/nEntries : 0.;
  //Event Loop
  for( int i = 0; i < nEntries; i++) {
    tree->GetEntry(i);if(i%500 == 0) cout << "Entry " << i << " of " << nEntries << endl;
    //aaa if(nL3 > 0) {
    if(nL3==0) continue;
    if(nL3 > (*l3AssociationMyBit).size() ) continue;
    for (int iMu = 0; iMu < nMu; iMu++) {
      if ( (*l3IsAssociated).at(iMu) == 0) continue;
      if( (*muAllGlobalMuons).at(iMu) == 0) continue;
      if(iMu>=nL3) continue;

      //cout << iMu << " of " << nMu << " and nL3 " << nL3 << " has MyBit " ; //<< (*l3AssociationMyBit).at(iMu) << endl;

      //if((*l3AssociationMyBit).at(iMu) & 1<<0) cout << "0" << endl;
      //if((*l3AssociationMyBit).at(iMu) & 1<<1) cout << "1" << endl;
      //if((*l3AssociationMyBit).at(iMu) & 1<<2) cout << "2" << endl;
      //if((*l3AssociationMyBit).at(iMu) & 1<<3) cout << "3" << endl;
      //if((*l3AssociationMyBit).at(iMu) & 1<<4) cout << "4" << endl;
      //if((*l3AssociationMyBit).at(iMu) & 1<<5) cout << "5" << endl;

      l3TestPtAll->Fill((*l3Pt).at(iMu));

      for(unsigned int w=0 ; w<8 ; w++){
	int ii = w;

      	if( (*l3AssociationMyBit).at(iMu) & 1<<ii ) {

	  mu_TMLastStationLoose[ii]->Fill((*muTMLastStationLoose).at(iMu),weight);
	  mu_TMLastStationTight[ii]->Fill((*muTMLastStationTight).at(iMu),weight);
	  mu_TM2DCompatibilityLoose[ii]->Fill((*muTM2DCompatibilityLoose).at(iMu),weight);
	  mu_TM2DCompatibilityTight[ii]->Fill((*muTM2DCompatibilityTight).at(iMu),weight);
	  mu_GlobalMuonPromptTight[ii]->Fill((*muGlobalMuonPromptTight).at(iMu),weight);
	  glb_pt[ii]->Fill((*l3Pt).at(iMu),weight);
	  glb_deltaPt[ii]->Fill((*l2Pt).at(iMu)-(*tkTrackPt).at(iMu),weight);
	  glb_calComp[ii]->Fill((*muCaloCompatibility).at(iMu),weight);
	  glb_segComp[ii]->Fill((*muSegmentCompatibility).at(iMu),weight);
	  glb_nChamber[ii]->Fill((*muNumberOfChambers).at(iMu),weight);
	  glb_nChamberMatch[ii]->Fill((*muNumberOfMatches).at(iMu),weight);
	  glb_muTrkKink[ii]->Fill((*muTrkKink).at(iMu),weight);
	  glb_muGlbKink[ii]->Fill((*muGlbKink).at(iMu),weight);
	  glb_muTrkRelChi2[ii]->Fill((*muTrkRelChi2).at(iMu),weight);
	  glb_muStaRelChi2[ii]->Fill((*muStaRelChi2).at(iMu),weight);

	  int dof = ((*tkTrackNdof).at(iMu)) ?  ((*tkTrackNdof).at(iMu)) : 1;
	  glb_muDeltaTrkRelChi2[ii]->Fill(fabs((*muTrkRelChi2).at(iMu) - (*tkTrackChi2).at(iMu)/dof ),weight);
	  dof = ((*l2Ndof).at(iMu)) ?  ((*l2Ndof).at(iMu)) : 1;
	  glb_muDeltaStaRelChi2[ii]->Fill(fabs((*muStaRelChi2).at(iMu) - (*l2Chi2).at(iMu)/dof ),weight);
	  if((*muAllGlobalMuons).at(iMu)==1 && (*l3Ndof).at(iMu)>0 ) glb_nChi2[ii]->Fill((*l3Chi2).at(iMu)/(*l3Ndof).at(iMu),weight);

	  ////for (int station = 0; station < 4; station++) {
	  double CSC1=CSC2=CSC3=CSC4=0.;
	  for (map<int,std::vector<int> >::const_iterator ggg = (*muNCSCSeg).begin(); ggg != (*muNCSCSeg).end(); ggg++) {        
	    CSC1 = ggg->second.at(0);
	    CSC2 = ggg->second.at(1);
	    CSC3 = ggg->second.at(2);
	    CSC4 = ggg->second.at(3);
	    glb_NCSC_1[ii]->Fill(ggg->second.at(0),weight);
	    glb_NCSC_2[ii]->Fill(ggg->second.at(1),weight);
	    glb_NCSC_3[ii]->Fill(ggg->second.at(2),weight);
	    glb_NCSC_4[ii]->Fill(ggg->second.at(3),weight);
	  }
	  glb_NCSC_total[ii]->Fill(CSC1+CSC2+CSC3+CSC4,weight);
	  double DT1=DT2=DT3=DT4=0.;
	  for (map<int,std::vector<int> >::const_iterator ggg = (*muNDTSeg).begin(); ggg != (*muNDTSeg).end(); ggg++) {        
	    DT1 = ggg->second.at(0);
	    DT2 = ggg->second.at(1);
	    DT3 = ggg->second.at(2);
	    DT4 = ggg->second.at(3);
	    glb_NDT_1[ii]->Fill(ggg->second.at(0),weight);
	    glb_NDT_2[ii]->Fill(ggg->second.at(1),weight);
	    glb_NDT_3[ii]->Fill(ggg->second.at(2),weight);
	    glb_NDT_4[ii]->Fill(ggg->second.at(3),weight);	    
	  }
	  glb_NDT_total[ii]->Fill(DT1+DT2+DT3+DT4,weight);
	  ////}
	  //for(map<int,std::vector<int> >::const_iterator hits = (*l3MuStationNumber)->begin(,weight); hits != (*l3MuStationNumber)->end(,weight); hits++){
	  //if(hits.first=iMu) cout << "iMu " << iMu << " hits last " << hits.second.last() << endl;
	  //}
	  //tom
	  /*
	  for (map<int,vector<double> >::const_iterator hitsX = (*l2RecHitsX).begin(); hitsX != (*l2RecHitsX).end(); hitsX++) {        
	    for (map<int,vector<double> >::const_iterator hitsY = (*l2RecHitsY).begin(); hitsY != (*l2RecHitsY).end(); hitsY++) {          
	      for (map<int,vector<double> >::const_iterator hitsZ = (*l2RecHitsZ).begin(); hitsZ != (*l2RecHitsZ).end(); hitsZ++) {
		if (hitsX->first == hitsY->first && hitsX->first == hitsZ->first) { //same L2 muon              
		  for (int i = 0; i < hitsX->second.size(); i++) {
		    poorMansDisplay->Fill(hitsX->second.at(i),hitsY->second.at(i),hitsZ->second.at(i));
		  }            
		}
	      } 
	    }
	  }
	  
	  for (map<int,std::vector<double> >::const_iterator hits = (*l2RecHitsX).begin(); hits != (*l2RecHitsX).end(); hits++) {        
	    for (int i = 0; i < hits->second.size(); i++) {          
	      recHitsTest->Fill(hits->second.at(i));       
	    } 
	  }
	  */	  
	//end tom
	}
      }
      if((*muAllGlobalMuons).at(iMu)==1) l3deltaPt->Fill((*l2Pt).at(iMu)-(*tkTrackPt).at(iMu),weight);
    }
    //aaa }
  }

  for(int j =0; j<8;j++){
    //h_pt[j]->Write("",TObject::kOverwrite);
    dirs[j]->Write("",TObject::kOverwrite);
  }

  histDirAlgo->Write("",TObject::kOverwrite);
  histFileAlgo->Close();

  return 0;
}
