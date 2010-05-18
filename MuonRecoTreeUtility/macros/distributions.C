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

  TH2F* mu_vtx[8];

  TH1F* glb_pt[8];
  TH1F* glb_deltaPtSTA[8];
  TH1F* glb_deltaPtTK[8];
  TH1F* glb_norm_deltaPt[8];
  TH1F* glb_norm_deltaPtSTA[8];
  TH1F* glb_norm_deltaPtTK[8];
  TH1F* glb_deltaPt[8];
  TH2F* glb_pt2D[8];
  TH1F* glb_calComp[8];
  TH1F* glb_segComp[8];
  TH1F* glb_nChamber[8];
  TH1F* glb_nChamberMatch[8];
  TH1F* glb_depth[8];
  TH1F* glb_NCSC_1[8];
  TH1F* glb_NCSC_2[8];
  TH1F* glb_NCSC_3[8];
  TH1F* glb_NCSC_4[8];
  TH1F* glb_NCSCArb_1[8];
  TH1F* glb_NCSCArb_2[8];
  TH1F* glb_NCSCArb_3[8];
  TH1F* glb_NCSCArb_4[8];
  TH1F* glb_NDT_1[8];
  TH1F* glb_NDT_2[8];
  TH1F* glb_NDT_3[8];
  TH1F* glb_NDT_4[8];
  TH1F* glb_NDTArb_1[8];
  TH1F* glb_NDTArb_2[8];
  TH1F* glb_NDTArb_3[8];
  TH1F* glb_NDTArb_4[8];
  TH1F* glb_NRPC_1[8];
  TH1F* glb_NRPC_2[8];
  TH1F* glb_NRPC_3[8];
  TH1F* glb_NRPC_4[8];
  TH1F* glb_NRPCArb_1[8];
  TH1F* glb_NRPCArb_2[8];
  TH1F* glb_NRPCArb_3[8];
  TH1F* glb_NRPCArb_4[8];
  TH1F* glb_NCSC_total[8];
  TH1F* glb_NCSCArb_total[8];
  TH1F* glb_NDT_total[8];
  TH1F* glb_NDTArb_total[8];
  TH1F* glb_NRPC_total[8];
  TH1F* glb_NRPCArb_total[8];
  TH1F* glb_NSeg_total[8];
  TH1F* glb_NSegArb_total[8];
  TH1F* glb_NHits[8];
  TH1F* glb_NMuHits[8];
  TH1F* glb_NDTHits[8];
  TH1F* glb_NCSCHits[8];
  TH1F* glb_NRPCHits[8];
  TH1F* glb_hitsIn1Ratio[8];
  TH2F* glb_hitsIn1RatioNMuHits[8];
  TH2F* glb_hitsIn1RatioDepth[8];
  TH2F* glb_hitsIn1RatioPt[8];
  TH1F* glb_muTrkKink[8];
  TH1F* glb_muGlbKink[8];
  TH1F* glb_muTrkRelChi2[8];
  TH1F* glb_muStaRelChi2[8];
  TH1F* glb_muDeltaTrkRelChi2[8];
  TH1F* glb_muDeltaStaRelChi2[8];
  TH1F* glb_nChi2[8];

  TH1F* glb_tkIso[8];
  TH1F* glb_tkIsoNTrk[8];
  TH1F* glb_normTkIso[8];
  TH1F* glb_calIso[8];

  TH1F* glb_d0[8];
  TH1F* sta_d0[8];
  TH1F* trk_d0[8];

  TH1F* mu_TMLastStationLoose[8];
  TH1F* mu_TMLastStationTight[8];
  TH1F* mu_TM2DCompatibilityLoose[8];
  TH1F* mu_TM2DCompatibilityTight[8];
  TH1F* mu_GlobalMuonPromptTight[8];
  TH1F* glb_minimumRequirement[8];
  TH1F* glb_testRatioBool[8];

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

    mu_vtx[ww] = new TH2F("mu_vtx","Muon Decay Vertex",2000,-1000.,1000.,1000,0.,1000.);

    glb_pt2D[ww] = new TH2F("glb_pt2D","histogram of p_{T}",100,0,500,100,0,500);
    glb_deltaPt[ww] = new TH1F("deltaPt","histogram of #Delta p_{T}",100,0,50);
    glb_deltaPtSTA[ww] = new TH1F("deltaPtSTA","histogram of #Delta p_{T}",100,0,50);
    glb_deltaPtTK[ww] = new TH1F("deltaPtTK","histogram of #Delta p_{T}",100,0,50);
    glb_norm_deltaPt[ww] = new TH1F("normDeltaPt","histogram of #Delta p_{T}",100,0,5);
    glb_norm_deltaPtSTA[ww] = new TH1F("normDeltaPtSTA","histogram of #Delta p_{T}",100,0,5);
    glb_norm_deltaPtTK[ww] = new TH1F("normDeltaPtTK","histogram of #Delta p_{T}",100,0,5);
    mu_TMLastStationLoose[ww] = new TH1F("mu_TMLastStationLoose","mu_TMLastStationLoose",3,-1.5,1.5);
    mu_TMLastStationTight[ww] = new TH1F("mu_TMLastStationTight","mu_TMLastStationTight",3,-1.5,1.5);;
    mu_TM2DCompatibilityLoose[ww] = new TH1F("mu_TM2DCompatibilityLoose","mu_TM2DCompatibilityLoose",3,-1.5,1.5);;
    mu_TM2DCompatibilityTight[ww] = new TH1F("mu_TM2DCompatibilityTight","mu_TM2DCompatibilityTight",3,-1.5,1.5);;
    mu_GlobalMuonPromptTight[ww] = new TH1F("mu_GlobalMuonPromptTight","mu_GlobalMuonPromptTight",3,-1.5,1.5);;
    glb_minimumRequirement[ww] = new TH1F("glb_minimumRequirement","Only 1 DT Segment and 1 RPC Hit",3,-1.5,1.5);;
    glb_testRatioBool[ww] = new TH1F("glb_testRatioBool","Test Ratio",3,-1.5,1.5);;
    glb_calComp[ww] = new TH1F("glb_calComp","Calo Compatibility",20,0.,1.);
    glb_segComp[ww] = new TH1F("glb_segComp","Segment Compatibility",20,0.,1.);
    glb_nChamber[ww] = new TH1F("glb_nChamber","Number of Chambers",11,-0.5,10.5);
    glb_nChamberMatch[ww] = new TH1F("glb_nChamberMatch","Number of Matched Chambers",11,-0.5,10.5);
    glb_depth[ww] = new TH1F("glb_depth","Station Depth",6,-0.5,5.5);
    glb_NCSC_1[ww] = new TH1F("glb_nCSC_1","Number of segments in station 1",11,-0.5,10.5);
    glb_NCSC_2[ww] = new TH1F("glb_nCSC_2","Number of segments in station 2",11,-0.5,10.5);
    glb_NCSC_3[ww] = new TH1F("glb_nCSC_3","Number of segments in station 3",11,-0.5,10.5);
    glb_NCSC_4[ww] = new TH1F("glb_nCSC_4","Number of segments in station 4",11,-0.5,10.5);
    glb_NCSCArb_1[ww] = new TH1F("glb_nCSCArb_1","Number of segments in station 1",11,-0.5,10.5);
    glb_NCSCArb_2[ww] = new TH1F("glb_nCSCArb_2","Number of segments in station 2",11,-0.5,10.5);
    glb_NCSCArb_3[ww] = new TH1F("glb_nCSCArb_3","Number of segments in station 3",11,-0.5,10.5);
    glb_NCSCArb_4[ww] = new TH1F("glb_nCSCArb_4","Number of segments in station 4",11,-0.5,10.5);
    glb_NDT_1[ww] = new TH1F("glb_nDT_1","Number of segments in DT station 1",11,-0.5,10.5);
    glb_NDT_2[ww] = new TH1F("glb_nDT_2","Number of segments in DT station 2",11,-0.5,10.5);
    glb_NDT_3[ww] = new TH1F("glb_nDT_3","Number of segments in DT station 3",11,-0.5,10.5);
    glb_NDT_4[ww] = new TH1F("glb_nDT_4","Number of segments in DT station 4",11,-0.5,10.5);
    glb_NDTArb_1[ww] = new TH1F("glb_nDTArb_1","Number of segments in DT station 1",11,-0.5,10.5);
    glb_NDTArb_2[ww] = new TH1F("glb_nDTArb_2","Number of segments in DT station 2",11,-0.5,10.5);
    glb_NDTArb_3[ww] = new TH1F("glb_nDTArb_3","Number of segments in DT station 3",11,-0.5,10.5);
    glb_NDTArb_4[ww] = new TH1F("glb_nDTArb_4","Number of segments in DT station 4",11,-0.5,10.5);
    glb_NRPC_1[ww] = new TH1F("glb_nRPC_1","Number of segments in RPC station 1",11,-0.5,10.5);
    glb_NRPC_2[ww] = new TH1F("glb_nRPC_2","Number of segments in RPC station 2",11,-0.5,10.5);
    glb_NRPC_3[ww] = new TH1F("glb_nRPC_3","Number of segments in RPC station 3",11,-0.5,10.5);
    glb_NRPC_4[ww] = new TH1F("glb_nRPC_4","Number of segments in RPC station 4",11,-0.5,10.5);
    glb_NRPCArb_1[ww] = new TH1F("glb_nRPCArb_1","Number of segments in RPC station 1",11,-0.5,10.5);
    glb_NRPCArb_2[ww] = new TH1F("glb_nRPCArb_2","Number of segments in RPC station 2",11,-0.5,10.5);
    glb_NRPCArb_3[ww] = new TH1F("glb_nRPCArb_3","Number of segments in RPC station 3",11,-0.5,10.5);
    glb_NRPCArb_4[ww] = new TH1F("glb_nRPCArb_4","Number of segments in RPC station 4",11,-0.5,10.5);
    glb_NCSC_total[ww] = new TH1F("glb_nCSC_total","Number of segments in all CSC stations",11,-0.5,10.5);
    glb_NCSCArb_total[ww] = new TH1F("glb_nCSCArb_total","Number of segments in all CSC stations",11,-0.5,10.5);
    glb_NDT_total[ww] = new TH1F("glb_nDT_total","Number of segments in all DT stations",11,-0.5,10.5);
    glb_NDTArb_total[ww] = new TH1F("glb_nDTArb_total","Number of segments in all DT stations",11,-0.5,10.5);
    glb_NRPC_total[ww] = new TH1F("glb_nRPC_total","Number of segments in all RPC stations",11,-0.5,10.5);
    glb_NRPCArb_total[ww] = new TH1F("glb_nRPCArb_total","Number of segments in all RPC stations",11,-0.5,10.5);
    glb_NSeg_total[ww] = new TH1F("glb_nSeg_total","Number of segments in all stations",16,-0.5,15.5);
    glb_NSegArb_total[ww] = new TH1F("glb_nSegArb_total","Number of segments in all stations",16,-0.5,15.5);
    glb_NMuHits[ww] = new TH1F("glb_nMuHits","Number of muon hits",76,-0.5,75.5);
    glb_NDTHits[ww] = new TH1F("glb_nDTHits","Number of DT hits",76,-0.5,75.5);
    glb_NCSCHits[ww] = new TH1F("glb_nCSCHits","Number of CSC hits",76,-0.5,75.5);
    glb_NRPCHits[ww] = new TH1F("glb_nRPCHits","Number of RPC hits",76,-0.5,75.5);
    glb_NHits[ww] = new TH1F("glb_nHits","Number of hits",76,-0.5,75.5);
    glb_hitsIn1Ratio[ww] = new TH1F("glb_hitsIn1Ratio","Ratio of Station1 Hits",100,0.,1.);
    glb_hitsIn1RatioNMuHits[ww] = new TH2F("glb_hitsIn1RatioNMuHits","Ratio of Station1 Hits",76,-0.5,75.5,100,0.,1.);
    glb_hitsIn1RatioDepth[ww] = new TH2F("glb_hitsIn1RatioDepth","Ratio of Station1 Hits",6,-0.5,5.5,100,0.,1.);
    glb_hitsIn1RatioPt[ww] = new TH2F("glb_hitsIn1RatioPt","Ratio of Station1 Hits",250,0.,500,100,0.,1.);
    glb_muTrkKink[ww] = new TH1F("glb_trkKink","Tracker Kink",500,0,500);
    glb_muGlbKink[ww] = new TH1F("glb_glbKink","Global Kink",500,0,500);
    glb_muTrkRelChi2[ww] = new TH1F("glb_trkRelChi2","Tracker Relative #Chi^{2}",50,0,50);
    glb_muStaRelChi2[ww] = new TH1F("glb_staRelChi2","Standalone Relative #Chi^{2}",50,0,50);
    glb_muDeltaTrkRelChi2[ww] = new TH1F("glb_trkDeltaRelChi2","#Delta Tracker Relative #Chi^{2}",10,0,5);
    glb_muDeltaStaRelChi2[ww] = new TH1F("glb_staDeltaRelChi2","#Delta Standalone Relative #Chi^{2}",10,0,5);
    glb_nChi2[ww] = new TH1F("glb_nChi2","Global #Chi^{2}/DoF",50,0,50);

    glb_tkIso[ww] = new TH1F("glb_tkIso","Tracker Isolation",100,0,20);
    glb_tkIsoNTrk[ww] = new TH1F("glb_tkIsoNTrk","Tracker Isolation N Tracks",110,-0.5,10.5);
    glb_normTkIso[ww] = new TH1F("glb_normTkIso","Normalized Tracker Isolation",100,0,20);
    glb_calIso[ww] = new TH1F("glb_calIso","Calo Isolation",100,0,20);

    glb_d0[ww] = new TH1F("glb_d0","Global impact parameter",200,0,2);
    sta_d0[ww] = new TH1F("sta_d0","Stand-alone impact parameter",200,0,2);
    trk_d0[ww] = new TH1F("trk_d0","Tracker impact parameter",200,0,2);
  }
  histDirAlgo->cd();

  TH1F *l3TestPtAll = new TH1F("l3TestPtAll","histogram of L3 p_{T}",250,0,500);
  TH1F *l3TestPtSel = new TH1F("l3TestPtSel","histogram of L3 p_{T}",250,0,500);
  TH1F *tkTestPtAll = new TH1F("tkTestPtAll","histogram of Tk p_{T}",250,0,500);
  TH1F *tkTestPtSel = new TH1F("tkTestPtSel","histogram of Tk p_{T}",250,0,500);
  TH1F *l3deltaPt = new TH1F("l3deltaPt","histogram of L3 #Delta p_{T}",100,-50.,50.);

  //TH1F *recHitsTest = new TH1F("recHitsTest","test to fill x position of recHits",1600,-800,800); 
  //TH3F *poorMansDisplay = new TH3F("poorMansDisplay","test for hits", 200, -800.0, 800.0, 200, -800.0, 800.0, 250, -1045.0, 1045.0); 

  int nEntries = tree->GetEntries();
  //double weight = nEntries > 0 ? effCS/nEntries : 0.;
  double weight = nEntries > 0 ? effCS : 0.;
  //Event Loop
  for( int i = 0; i < nEntries; i++) {
    tree->GetEntry(i);if(i%500 == 0) cout << "Entry " << i << " of " << nEntries << endl;
    //ccc aaa if(nL3 > 0) {
    //aaa if(nL3==0) continue;
    //aaa if(nSimMuon>0) continue;
    //aaa if(nL3 > (*l3AssociationMyBit).size() ) continue;
    for (int iMu = 0; iMu < nMu; iMu++) {
      l3TestPtAll->Fill((*l3Pt).at(iMu));
      tkTestPtAll->Fill((*tkTrackPt).at(iMu));
      //aaa if ( (*l3IsAssociated).at(iMu) == 0) continue;
      //aaa if(nSimMuon>0) continue;
      //aaa if( (*muAllGlobalMuons).at(iMu) == 1) continue;
      //aaa if(iMu>=nL3) continue;
      //aaa if(fabs((*l3Eta).at(iMu))>0.8) continue;
      l3TestPtSel->Fill((*l3Pt).at(iMu));
      tkTestPtSel->Fill((*tkTrackPt).at(iMu));

      //adam temp depth
      unsigned int mask = (*muStationMask).at(iMu);
      bool station1 = ((mask & 1<<0)||(mask & 1<<4));
      bool station2 = ((mask & 1<<1)||(mask & 1<<5));
      bool station3 = ((mask & 1<<2)||(mask & 1<<6));
      bool station4 = ((mask & 1<<3)||(mask & 1<<7));
      
      bool depth1 = station1 && !station2 && !station3 && !station4;
      bool depth2 = station2 && !station3 && !station4;
      bool depth3 = station3 && !station4;
      bool depth4 = station4;
      //if(depth2 || depth3 || depth4) continue;
      //aaa if(fabs((*l3Eta).at(iMu))>0.8) continue;
      //end adam temp depth

      //cout << iMu << " of " << nMu << " and nL3 " << nL3 << " has MyBit " ; //<< (*l3AssociationMyBit).at(iMu) << endl;

      //if((*l3AssociationMyBit).at(iMu) & 1<<0) cout << "0" << endl;
      //if((*l3AssociationMyBit).at(iMu) & 1<<1) cout << "1" << endl;
      //if((*l3AssociationMyBit).at(iMu) & 1<<2) cout << "2" << endl;
      //if((*l3AssociationMyBit).at(iMu) & 1<<3) cout << "3" << endl;
      //if((*l3AssociationMyBit).at(iMu) & 1<<4) cout << "4" << endl;
      //if((*l3AssociationMyBit).at(iMu) & 1<<5) cout << "5" << endl;

      

      for(unsigned int w=0 ; w<8 ; w++){
	int ii = w;

      	if( (*l3AssociationMyBit).at(iMu) & 1<<ii ) {

	  float r = sqrt((*l3AssociationVtxX).at(iMu) 
			 * (*l3AssociationVtxX).at(iMu)
			 + (*l3AssociationVtxY).at(iMu) 
			 * (*l3AssociationVtxY).at(iMu));
	  float z = (*l3AssociationVtxZ).at(iMu);
	  mu_vtx[ii]->Fill(z,r,weight);

	  mu_TMLastStationLoose[ii]->Fill((*muTMLastStationLoose).at(iMu),weight);
	  mu_TMLastStationTight[ii]->Fill((*muTMLastStationTight).at(iMu),weight);
	  mu_TM2DCompatibilityLoose[ii]->Fill((*muTM2DCompatibilityLoose).at(iMu),weight);
	  mu_TM2DCompatibilityTight[ii]->Fill((*muTM2DCompatibilityTight).at(iMu),weight);
	  mu_GlobalMuonPromptTight[ii]->Fill((*muGlobalMuonPromptTight).at(iMu),weight);
	  glb_pt[ii]->Fill((*l3Pt).at(iMu),weight);
	  glb_pt2D[ii]->Fill((*tkTrackPt).at(iMu),(*l2Pt).at(iMu),weight);

	  if((*muIso03Valid).at(iMu)>0) glb_tkIso[ii]->Fill(((*muIso03sumPt).at(iMu)-(*muIso03trackerVetoPt).at(iMu)),weight);
	  if((*muIso03Valid).at(iMu)>0) glb_normTkIso[ii]->Fill((((*muIso03sumPt).at(iMu)-(*muIso03trackerVetoPt).at(iMu))/(*l3Pt).at(iMu)),weight);
	  if((*muIso03Valid).at(iMu)>0) glb_calIso[ii]->Fill(((*muIso03emEt).at(iMu)+(*muIso03hadEt).at(iMu)), weight);
	  if((*muIso03Valid).at(iMu)>0) glb_tkIsoNTrk[ii]->Fill((*muIso03nTracks).at(iMu),weight);
	  
	  glb_d0[ii]->Fill(-1.*(*l3D0).at(iMu),weight);
	  sta_d0[ii]->Fill(-1.*(*l2D0).at(iMu),weight);
	  trk_d0[ii]->Fill(-1.*(*tkTrackD0).at(iMu),weight);

	  glb_deltaPt[ii]->Fill(fabs((*tkTrackPt).at(iMu)-(*l2Pt).at(iMu)),weight);
	  glb_norm_deltaPt[ii]->Fill(fabs(((*tkTrackPt).at(iMu)-(*l2Pt).at(iMu))/(*tkTrackPt).at(iMu)),weight);

	  if((*l2Pt).at(iMu) >= (*tkTrackPt).at(iMu))
	    glb_deltaPtSTA[ii]->Fill(fabs((*l2Pt).at(iMu)-(*tkTrackPt).at(iMu)),weight);
	  if((*l2Pt).at(iMu) < (*tkTrackPt).at(iMu))
	    glb_deltaPtTK[ii]->Fill(fabs((*tkTrackPt).at(iMu)-(*l2Pt).at(iMu)),weight);
	  if((*l2Pt).at(iMu) >= (*tkTrackPt).at(iMu))
	    glb_norm_deltaPtSTA[ii]->Fill(fabs(((*l2Pt).at(iMu)-(*tkTrackPt).at(iMu)) / (*l2Pt).at(iMu)),weight);
	  if((*l2Pt).at(iMu) < (*tkTrackPt).at(iMu))
	    glb_norm_deltaPtTK[ii]->Fill(fabs(((*l2Pt).at(iMu)-(*tkTrackPt).at(iMu)) / (*tkTrackPt).at(iMu)),weight);

	  glb_calComp[ii]->Fill((*muCaloCompatibility).at(iMu),weight);
	  glb_segComp[ii]->Fill((*muSegmentCompatibility).at(iMu),weight);
	  glb_nChamber[ii]->Fill((*muNumberOfChambers).at(iMu),weight);
	  if( (*muNumberOfChambers).at(iMu) > 1 /*(*muTMOneStationLoose).at(iMu) > 0*/ ) glb_nChamberMatch[ii]->Fill((*muNumberOfMatches).at(iMu),weight);

	  unsigned int mask = (*muStationMask).at(iMu);

	  bool station1 = ((mask & 1<<0)||(mask & 1<<4));
	  bool station2 = ((mask & 1<<1)||(mask & 1<<5));
	  bool station3 = ((mask & 1<<2)||(mask & 1<<6));
	  bool station4 = ((mask & 1<<3)||(mask & 1<<7));
	  
	  bool depth1 = station1 && !station2 && !station3 && !station4;
	  bool depth2 = station2 && !station3 && !station4;
	  bool depth3 = station3 && !station4;
	  bool depth4 = station4;

	  if(depth1) glb_depth[ii]->Fill(1,weight);
	  if(depth2) glb_depth[ii]->Fill(2,weight);
	  if(depth3) glb_depth[ii]->Fill(3,weight);
	  if(depth4) glb_depth[ii]->Fill(4,weight);

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
	  double CSC1,CSC2,CSC3,CSC4;
	  CSC1 = CSC2 = CSC3 = CSC4 = 0;
	  for (map<int,std::vector<int> >::const_iterator ggg = (*muNCSCSeg).begin(); ggg != (*muNCSCSeg).end(); ggg++) {
	    if(ggg->first != iMu) continue;
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
	  double CSC1A,CSC2A,CSC3A,CSC4A;
	  CSC1A = CSC2A = CSC3A = CSC4A = 0;
	  for (map<int,std::vector<int> >::const_iterator ggg = (*muNCSCSegArb).begin(); ggg != (*muNCSCSegArb).end(); ggg++) {        
	    if (ggg->first != iMu) continue;
	    CSC1A = ggg->second.at(0);
	    CSC2A = ggg->second.at(1);
	    CSC3A = ggg->second.at(2);
	    CSC4A = ggg->second.at(3);
	    glb_NCSCArb_1[ii]->Fill(ggg->second.at(0),weight);
	    glb_NCSCArb_2[ii]->Fill(ggg->second.at(1),weight);
	    glb_NCSCArb_3[ii]->Fill(ggg->second.at(2),weight);
	    glb_NCSCArb_4[ii]->Fill(ggg->second.at(3),weight);
	  }
	  glb_NCSCArb_total[ii]->Fill(CSC1A+CSC2A+CSC3A+CSC4A,weight);
	  double DT1,DT2,DT3,DT4;
	  DT1 = DT2 = DT3 = DT4 = 0;
	  for (map<int,std::vector<int> >::const_iterator ggg = (*muNDTSeg).begin(); ggg != (*muNDTSeg).end(); ggg++) {        
	    if(ggg->first != iMu) continue;
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
	  double DT1A,DT2A,DT3A,DT4A;
	  DT1A = DT2A = DT3A = DT4A = 0;
	  for (map<int,std::vector<int> >::const_iterator ggg = (*muNDTSegArb).begin(); ggg != (*muNDTSegArb).end(); ggg++) {        
	    if(ggg->first != iMu) continue;
	    DT1A = ggg->second.at(0);
	    DT2A = ggg->second.at(1);
	    DT3A = ggg->second.at(2);
	    DT4A = ggg->second.at(3);
	    glb_NDTArb_1[ii]->Fill(ggg->second.at(0),weight);
	    glb_NDTArb_2[ii]->Fill(ggg->second.at(1),weight);
	    glb_NDTArb_3[ii]->Fill(ggg->second.at(2),weight);
	    glb_NDTArb_4[ii]->Fill(ggg->second.at(3),weight);	    
	  }
	  glb_NDTArb_total[ii]->Fill(DT1A+DT2A+DT3A+DT4A,weight);
	  double RPC1,RPC2,RPC3,RPC4;
	  RPC1 = RPC2 = RPC3 = RPC4 = 0;
	  for (map<int,std::vector<int> >::const_iterator ggg = (*muNRPCSeg).begin(); ggg != (*muNRPCSeg).end(); ggg++) {        
	    if(ggg->first != iMu) continue;
	    RPC1 = ggg->second.at(0);
	    RPC2 = ggg->second.at(1);
	    RPC3 = ggg->second.at(2);
	    RPC4 = ggg->second.at(3);
	    glb_NRPC_1[ii]->Fill(ggg->second.at(0),weight);
	    glb_NRPC_2[ii]->Fill(ggg->second.at(1),weight);
	    glb_NRPC_3[ii]->Fill(ggg->second.at(2),weight);
	    glb_NRPC_4[ii]->Fill(ggg->second.at(3),weight);	    
	  }
	  glb_NRPC_total[ii]->Fill(RPC1+RPC2+RPC3+RPC4,weight);
	  double RPC1A,RPC2A,RPC3A,RPC4A;
	  RPC1A = RPC2A = RPC3A = RPC4A = 0;
	  for (map<int,std::vector<int> >::const_iterator ggg = (*muNRPCSegArb).begin(); ggg != (*muNRPCSegArb).end(); ggg++) {        
	    if(ggg->first != iMu) continue;
	    RPC1A = ggg->second.at(0);
	    RPC2A = ggg->second.at(1);
	    RPC3A = ggg->second.at(2);
	    RPC4A = ggg->second.at(3);
	    glb_NRPCArb_1[ii]->Fill(ggg->second.at(0),weight);
	    glb_NRPCArb_2[ii]->Fill(ggg->second.at(1),weight);
	    glb_NRPCArb_3[ii]->Fill(ggg->second.at(2),weight);
	    glb_NRPCArb_4[ii]->Fill(ggg->second.at(3),weight);	    
	  }
	  glb_NRPCArb_total[ii]->Fill(RPC1A+RPC2A+RPC3A+RPC4A,weight);
	  glb_NSeg_total[ii]->Fill(CSC1+CSC2+CSC3+CSC4+DT1+DT2+DT3+DT4,weight);
	  glb_NSegArb_total[ii]->Fill(CSC1A+CSC2A+CSC3A+CSC4A+DT1A+DT2A+DT3A+DT4A,weight);

	  //int nMuHits = 0;
	  glb_NHits[ii]->Fill((*l3NHits).at(iMu),weight);
	  glb_NMuHits[ii]->Fill((*l3NMuHits).at(iMu),weight);
	  glb_NDTHits[ii]->Fill((*l3NDTHits).at(iMu),weight);
	  glb_NCSCHits[ii]->Fill((*l3NCSCHits).at(iMu),weight);
	  glb_NRPCHits[ii]->Fill((*l3NRPCHits).at(iMu),weight);
	  /*
	  for (map<int, int>::const_iterator hits = (*l3NMuHits).begin(); hits != (*l3NMuHits).end(); hits++) {
	    if(hits->first==iMu && hits->second > 0) glb_NMuHits[ii]->Fill(hits->second,weight);
	    //nMuHits = hits->second;
	  }
	  for (map<int, int>::const_iterator hits = (*l3NDTHits).begin(); hits != (*l3NDTHits).end(); hits++) {
	    if(hits->first==iMu && hits->second > 0) glb_NDTHits[ii]->Fill(hits->second,weight);
	    //nMuHits = hits->second;
	  }
	  for (map<int, int>::const_iterator hits = (*l3NCSCHits).begin(); hits != (*l3NCSCHits).end(); hits++) {
	    if(hits->first==iMu && hits->second > 0) glb_NCSCHits[ii]->Fill(hits->second,weight);
	    //nMuHits = hits->second;
	  }
	  for (map<int, int>::const_iterator hits = (*l3NRPCHits).begin(); hits != (*l3NRPCHits).end(); hits++) {
	    if(hits->first==iMu && hits->second > 0) glb_NRPCHits[ii]->Fill(hits->second,weight);
	    //nMuHits = hits->second;
	  }
	  */
	  
	  float hitsIn1 = 0.;
	  float hitsIn2 = 0.;
	  float hitsIn3 = 0.;
	  float hitsIn4 = 0.;
	  float hitsInAll = 0.;
	  for (map<int,std::vector<int> >::const_iterator hits1 = (*l3MuStationNumber).begin(); hits1 != (*l3MuStationNumber).end(); hits1++) {
	    if(hits1->first == iMu) {
	      for(int i1 = 0; i1 < hits1->second.size(); i1++) {
		if(hits1->second.at(i1) == 1) hitsIn1 += 1.0;
		if(hits1->second.at(i1) == 2) hitsIn2 += 1.0;
		if(hits1->second.at(i1) == 3) hitsIn3 += 1.0;
		if(hits1->second.at(i1) == 4) hitsIn4 += 1.0;
		hitsInAll += 1.0;
		// cout << "event " << i << " iMu " << iMu << " hits " << i1 << " " << hitsIn1 << " " << hitsInAll << endl;
	      }
	    }
	  }
	  if(hitsInAll>0 && hitsIn1>0) {
	    float ratio = hitsIn1/hitsInAll ;
	    glb_hitsIn1Ratio[ii]->Fill(ratio,weight);
	    glb_hitsIn1RatioNMuHits[ii]->Fill(hitsInAll,ratio,weight);
	    if(depth1) glb_hitsIn1RatioDepth[ii]->Fill(1,ratio,weight);
	    if(depth2) glb_hitsIn1RatioDepth[ii]->Fill(2,ratio,weight);
	    if(depth3) glb_hitsIn1RatioDepth[ii]->Fill(3,ratio,weight);
	    if(depth4) glb_hitsIn1RatioDepth[ii]->Fill(4,ratio,weight);
	    glb_hitsIn1RatioPt[ii]->Fill((*l3Pt).at(iMu),ratio,weight);
	    if(hitsInAll>10 && ratio < 0.94) {
	      glb_testRatioBool[ii]->Fill(1,weight);
	    } else {
	      glb_testRatioBool[ii]->Fill(0,weight);
	    }
	  } 

	  int rpcIn1 = 0;
	  for (map<int,std::vector<int> >::const_iterator hits1 = (*l3SubdetIds).begin(); hits1 != (*l3SubdetIds).end(); hits1++) {
	    if(hits1->first == iMu) {
	      for(int i1 = 0; i1 < hits1->second.size(); i1++) {
		if(hits1->second.at(i1) == 3) rpcIn1++;
	      }
	    }
	  }
	  if(DT1+DT2+DT3+DT4>0) {
	    if(DT1==1 && DT2==0 && DT3==0 && DT4==0) {
	      glb_minimumRequirement[ii]->Fill(1,weight);
	    } else {
	      glb_minimumRequirement[ii]->Fill(0,weight);
	    }
	  }
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
