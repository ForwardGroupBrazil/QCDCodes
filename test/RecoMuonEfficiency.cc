// -*- C++ -*-
//
// Package:    RecoMuonEfficiency
// Class:      RecoMuonEfficiency
// 
/**\class RecoMuonEfficiency RecoMuonEfficiency.cc test/RecoMuonEfficiency/src/RecoMuonEfficiency.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam Everett
//         Created:  Tue Sep  5 13:45:40 CDT 2006
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

//#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "CLHEP/HepMC/ReadHepMC.h"
#include "CLHEP/HepMC/GenEvent.h"

// algorithm include files
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TMath.h>
#include <TF1.h>

//
// class decleration
//

using namespace edm;
using namespace std;

class RecoMuonEfficiency : public edm::EDAnalyzer {
   public:
      explicit RecoMuonEfficiency(const edm::ParameterSet&);
      ~RecoMuonEfficiency();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual TH1F* divideErr(TH1F*, TH1F*);
      virtual float ptL2_90(float pt, float sigma_pt);
      virtual float ptL3_90(float pt, float sigma_pt);
      // ----------member data ---------------------------
  //used to select what tracks to read from configuration file
  edm::InputTag GLBtrackTags_; 
  edm::InputTag STAtrackTags_; 
  edm::InputTag TKtrackTags_; 
  edm::InputTag G4trackTags_; 
  edm::InputTag simVtxTags_; 

  MuonServiceProxy* theService;

  string out, open;
  double  etamin, etamax, ptmin;
  int nbins, partId;

  //const int nbins, partID;
  //const float etamin;
  //const float etamax;
  //const float ptmin ;
  
  bool online;

  //ROOT Pointers
  TFile* hFile;
  TStyle* effStyle;
  // efficiency 
  TH1F* hi_ptgen  ;
  TH1F* hi_ptL1   ;
  TH1F* hi_ptL2   ;
  TH1F* hi_ptL3   ;
  TH1F* hi_ptL2a  ;
  TH1F* hi_ptL3a  ;
     
  TH1F* hi_etagen ;
  TH1F* hi_etaL1  ;
  TH1F* hi_etaL2  ;
  TH1F* hi_etaL3  ;
  TH1F* hi_etaL2a ;
  TH1F* hi_etaL3a ;

  TH2F* hi_effgen ;
  TH2F* hi_effL1  ;
  TH2F* hi_effL2  ;
  TH2F* hi_effL3  ;

  TH1F* hi_ra ;
  TH1F* hi_l1 ;
  TH1F* hi_l1a ;
  TH1F* hi_l2 ;
  TH1F* hi_l3 ;

  TH1F* hi_l110pt ;
  TH1F* hi_l120pt ;
  TH1F* hi_l130pt ;
  TH1F* hi_l140pt ;

  TH1F* hi_l210pt ;
  TH1F* hi_l220pt ;
  TH1F* hi_l230pt ;
  TH1F* hi_l240pt ;

  TH1F* hi_l310pt ;
  TH1F* hi_l320pt ;
  TH1F* hi_l330pt ;
  TH1F* hi_l340pt ;

  // pt resolution
  TH1F* hi_ptres_L1  ;
  TH1F* hi_ptres_L2  ;
  TH1F* hi_ptres_L3  ;
  
  TH1F* hi_ptres_L1b ;
  TH1F* hi_ptres_L2b ;
  TH1F* hi_ptres_L3b ;
  TH1F* hi_ptres_L1o ;
  TH1F* hi_ptres_L2o ;
  TH1F* hi_ptres_L3o ;
  TH1F* hi_ptres_L1e ;
  TH1F* hi_ptres_L2e ;
  TH1F* hi_ptres_L3e ;

  // chi2
  TH1F* hi_chi2_L2 ;
  TH1F* hi_chi2_L3 ;
  TH1F* hi_prob_L2 ;
  TH1F* hi_prob_L3 ;
  
  // L2 pulls
  TH1F* hi_ptpull_L2 ;
  TH1F* hi_thetapull_L2 ;
  TH1F* hi_phipull_L2 ;

  // L3 pulls
  TH1F* hi_ptpull_L3 ;
  TH1F* hi_thetapull_L3 ;
  TH1F* hi_phipull_L3 ;

  // kink
  TH1F* hi_kink   ;

  TH1F* hi_regionR;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RecoMuonEfficiency::RecoMuonEfficiency(const edm::ParameterSet& pset)
:
  GLBtrackTags_(pset.getUntrackedParameter<edm::InputTag>("GLBtracks")),
  STAtrackTags_(pset.getUntrackedParameter<edm::InputTag>("STAtracks")),
  TKtrackTags_(pset.getUntrackedParameter<edm::InputTag>("TKtracks")),
  G4trackTags_(pset.getUntrackedParameter<edm::InputTag>("G4tracks")),
  simVtxTags_(pset.getUntrackedParameter<edm::InputTag>("simVtx")),
  out(pset.getParameter<string>("out")),
  open(pset.getParameter<string>("open")),
  etamin(pset.getParameter<double>("etamin")),
  etamax(pset.getParameter<double>("etamax")),
  ptmin(pset.getParameter<double>("ptmin")),
  nbins(pset.getParameter<int>("nbins")),
  partId(pset.getParameter<int>("partId"))
{
   //now do what ever initialization is needed
 
  // service parameters
  ParameterSet serviceParameters = pset.getParameter<ParameterSet>("ServiceParameters");
  // the services
  theService = new MuonServiceProxy(serviceParameters);
  
    hFile = new TFile( out.c_str(), open.c_str() );
}


RecoMuonEfficiency::~RecoMuonEfficiency()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  if (hFile!=0) {
    hFile->Close();
    delete hFile;
  }
  if (theService) delete theService;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
RecoMuonEfficiency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using reco::TrackCollection;

   // Update the services
   theService->update(iSetup);
  
   Handle<TrackCollection> GLBTrackCollection;
   iEvent.getByLabel( GLBtrackTags_, GLBTrackCollection);
   const reco::TrackCollection glbTC = *(GLBTrackCollection.product());

   Handle<TrackCollection> STATrackCollection;
   iEvent.getByLabel( STAtrackTags_, STATrackCollection);
   const reco::TrackCollection staTC = *(STATrackCollection.product());

   Handle<TrackCollection> TKTrackCollection;
   iEvent.getByLabel(TKtrackTags_,TKTrackCollection);
   const reco::TrackCollection tkTC = *(TKTrackCollection.product());   

   cout << "Tk: " << tkTC.size() << " STA: " << staTC.size() << " GLB: " << glbTC.size() << endl;

   //-----
   for (  TrackCollection::const_iterator staCand = staTC.begin(); staCand != staTC.end(); staCand++) {
     //TrajectoryStateOnSurface innerMuTsos;  
     //TrajectoryStateOnSurface outerTkTsos;
     
     //reco::TransientTrack staT(*staCand,&*theService->magneticField(),theService->trackingGeometry());  
     //innerMuTsos = staT.innermostMeasurementState();
     
     //float deltaR = 100.0;
     
     //int iPosition = 0;
     reco::TrackCollection::const_iterator is;
     for ( is = tkTC.begin(); is != tkTC.end(); ++is ) {
       //iPosition++;
       //double deltaEta = innerMuTsos.globalMomentum().eta() - is->eta();
       //double deltaPhi = innerMuTsos.globalMomentum().phi() - is->phi();
       double deltaEta = staCand->eta() - is->eta();
       double deltaPhi = staCand->phi() - is->phi();
       double deltaR_tmp = sqrt(pow(deltaEta,2.) + pow(deltaPhi,2.));
       
       hi_regionR->Fill(fabs(deltaR_tmp));
       
     }  
   }
   //-----
   
   Handle<SimTrackContainer> simTrackCollection;
   iEvent.getByLabel(G4trackTags_, simTrackCollection);
   const SimTrackContainer simTC = *(simTrackCollection.product());

   edm::Handle<SimVertexContainer> simVertexCollection;
   iEvent.getByLabel(simVtxTags_, simVertexCollection);
   const SimVertexContainer simVC = *(simVertexCollection.product());

   //edm::Handle< edm::HepMCProduct > hepmc ;
   //iEvent.getByLabel( "source", hepmc ) ;
   //const edm::HepMCProduct hepMC = *(hepmc.product());

   //const HepMC::GenEvent *hepmcEvent = hepMC.GetEvent();
   //for (HepMC::GenEvent::particle_iterator iter = hepmcEvent.particles_begin(); iter != hepmcEvent.particles_end(); ++iter) {
   //
   //}

   // weight
   const float LumiLHC = 1.;
   bool valid = false;
   int idx_mu = -1;
   float q = 0.0;
   float pt = 0.0;
   float eta = 0.0;
   float phi = 0.0;
   float vtx2d = 0.0;
   float vtx3d = 0.0;
   int orig = 0;    
   int region = -1;  // ( 0 = brl, 1 = ovl, 2 = fwd )
   
   float anaW = 1.0;
   
   //FIXME: this block should be GEN not SIM!!
   int l = 0;
   for (SimTrackContainer::const_iterator simTrack = simTC.begin(); simTrack != simTC.end(); simTrack++) {
     float eta1 = fabs(simTrack->momentum().pseudoRapidity());
     float pt1 = simTrack->momentum().perp();
     if ( eta1 >= etamin && eta1 <= etamax && pt1 >= ptmin ) {
       if ( pt1 >= pt ) {
	 idx_mu = l;
	 pt = pt1;
	 eta = simTrack->momentum().pseudoRapidity();
	 phi = simTrack->momentum().phi();
	 //FIXME : simtrack->charge()?
	 //q = simTrack->charge();//t1.Chagen[l];
	 q = -1;

	 if ( fabs(eta) <= 0.8 ) region = 0;
	 if ( fabs(eta) > 0.8 && fabs(eta) < 1.2 ) region = 1;
	 if ( fabs(eta) >= 1.2 ) region = 2;
       }
     }
     l++;
   }
   
   if ( pt > ptmin ) valid = true;
   
   if ( valid ) {
     
     bool L1valid = false;
     bool L2valid = false;
     bool L3valid = false;
     
     //
     // L1 trigger       
     //
     int idx_L1 = -1;
     float ptL1 = 0.0;
     float etaL1 = 0.0;
     float phiL1 = 0.0;
     float qL1 = 0.0;
     l=0;
     for (SimTrackContainer::const_iterator simTrack = simTC.begin(); simTrack != simTC.end(); simTrack++) {
       //ADAM: BxL1??    
       //if ( t1.BxL1[l] == 0 && t1.PtL1[l] >= ptL1 ) {
       if ( simTrack->momentum().perp()  >= ptL1 ) {
	 float eta1 = fabs(simTrack->momentum().pseudoRapidity()); 
	 if ( eta1 > etamax+0.02 || eta1 < etamin-0.02 ) continue;
	 if ( simTrack->momentum().perp() > ptL1 ) {
	   idx_L1 = l;
	   ptL1 = simTrack->momentum().perp();
	   etaL1 = simTrack->momentum().pseudoRapidity();
	   phiL1 = simTrack->momentum().phi();
	   //FIXME: simTrack->Charge()
	   qL1 = simTrack->type()/abs(simTrack->type());
	   //cout << "simCharge " << qL1 << endl;
	   //qL1 = simTrack->charge();//t1.ChaL1[l];
	   //qL1 = -1;
	 }
       }
       l++;
     }
     if ( idx_L1 > -1 ) L1valid = true;
     
     //
     // L2/standalone muon reconstruction
     //
     int idx_L2 = -1;
     float ptL2 = 0.0;
     float EptL2 = 0.0;
     float etaL2 = 0.0;
     float phiL2 = 0.0;
     float qL2 = 0.0;
     
     float ptvL2 = 0.0;
     float chi2L2 = 0.0;
     float probL2 = 0.0;
     
     {
       l=0;
       for (  TrackCollection::const_iterator staTrack = staTC.begin(); staTrack != staTC.end(); staTrack++) {
	 reco::TransientTrack theTT(*staTrack,&*theService->magneticField(),theService->trackingGeometry());
	 TrajectoryStateOnSurface firstState = theTT.impactPointState();
	 float ptv = firstState.globalMomentum().perp();
	 if ( ptv == 0.0 ) continue;
	 if ( ptv > ptL2 ) {
	   ptL2 = ptv;
	   //FIXME: EPtv?
	   //EptL2 = t1.EPtvSA[l];
	   etaL2 = firstState.globalMomentum().eta();
	   phiL2 = firstState.globalMomentum().phi();
	   qL2   = staTrack->charge();
	   cout << "staCharge " << qL2 << endl;
	   ptvL2 = ptv;
	   chi2L2 = (staTrack->ndof() > 0 ) ? staTrack->chi2()/staTrack->ndof() : 0;
	   probL2 =  TMath::Prob(staTrack->chi2(),(int)staTrack->ndof());
	   idx_L2 = l;
	 }
	 if ( idx_L2 > -1  ) L2valid = true;
       }
       l++;      
     }
     
     ptL2 = ptL2_90(ptL2,EptL2);
     
     //  
     // L3/global muon reconstruction
     //
     int idx_L3 = -1;
     float ptL3 = 0.0;
     float EptL3 = 0.0;
     float etaL3 = 0.0;
     float phiL3 = 0.0;
     float qL3 = 0.0;
     
     float ptvL3 = 0.0;
     float chi2L3 = 0.0;
     float probL3 = 0.0;
     float kink = 0.0;
     
     { 
       l=0; 
       for ( TrackCollection::const_iterator glbTrack = glbTC.begin(); glbTrack != glbTC.end(); glbTrack++) {
	 if ( glbTrack->pt() > ptL3 ) {
	   ptL3 = glbTrack->pt();
	   //FIXME: EPtv?
	   //EptL3 = glbTrack->pt();//t1.EPtvGL[l];
	   etaL3 = glbTrack->eta();
	   phiL3 = glbTrack->phi();
	   qL3 = glbTrack->charge();
	   ptvL3 = ptL3;
	   chi2L3 = glbTrack->chi2()/glbTrack->ndof();
	   probL3 =  TMath::Prob(glbTrack->chi2(),(int)glbTrack->ndof());
	   idx_L3 = l; 
	 }
	 if ( L2valid && idx_L3 > -1  ) L3valid = true;
       }
       l++;
     }
     
     ptL3 = ptL3_90(ptL3,EptL3);
     
     //
     // fill histograms
     //
     hi_ptgen->Fill(pt,anaW);
     hi_etagen->Fill(fabs(eta),anaW);
     hi_effgen->Fill(pt,eta,anaW);
     
     if ( L1valid ) {
       hi_ptL1->Fill(pt,anaW);
       hi_etaL1->Fill(fabs(eta),anaW);
       hi_effL1->Fill(pt,eta,anaW);
       hi_ptres_L1->Fill((qL1/ptL1-q/pt)/(q/pt),anaW);
       if ( region == 0 ) hi_ptres_L1b->Fill((qL1/ptL1-q/pt)/(q/pt),anaW);
       if ( region == 1 ) hi_ptres_L1o->Fill((qL1/ptL1-q/pt)/(q/pt),anaW);
       if ( region == 2 ) hi_ptres_L1e->Fill((qL1/ptL1-q/pt)/(q/pt),anaW);
     }
     
     if ( L1valid && ptL1 >= 10 ) hi_l110pt->Fill(pt,anaW/LumiLHC);
     if ( L1valid && ptL1 >= 20 ) hi_l120pt->Fill(pt,anaW/LumiLHC);
     if ( L1valid && ptL1 >= 30 ) hi_l130pt->Fill(pt,anaW/LumiLHC);
     if ( L1valid && ptL1 >= 40 ) hi_l140pt->Fill(pt,anaW/LumiLHC);
     
     if ( L2valid && ptL2 >= 10 && ptL1 >= 10 ) hi_l210pt->Fill(pt,anaW/LumiLHC);
     if ( L3valid && ptL3 >= 10 && ptL2 >= 10 && ptL1 >= 10 ) hi_l310pt->Fill(pt,anaW/LumiLHC);
     if ( L2valid && ptL2 >= 20 && ptL1 >= 20 ) hi_l220pt->Fill(pt,anaW/LumiLHC);
     if ( L3valid && ptL3 >= 20 && ptL2 >= 20 && ptL1 >= 20 ) hi_l320pt->Fill(pt,anaW/LumiLHC);
     if ( L2valid && ptL2 >= 30 && ptL1 >= 30 ) hi_l230pt->Fill(pt,anaW/LumiLHC);
     if ( L3valid && ptL3 >= 30 && ptL2 >= 30 && ptL1 >= 30 ) hi_l330pt->Fill(pt,anaW/LumiLHC);
     if ( L2valid && ptL2 >= 40 && ptL1 >= 40 ) hi_l240pt->Fill(pt,anaW/LumiLHC);
     if ( L3valid && ptL3 >= 40 && ptL2 >= 40 && ptL1 >= 40 ) hi_l340pt->Fill(pt,anaW/LumiLHC);
     
     if ( L2valid ) {
       hi_ptL2->Fill(pt,anaW);
       hi_etaL2->Fill(fabs(eta),anaW);
       hi_effL2->Fill(pt,eta,anaW);
       hi_ptL2a->Fill(pt,anaW);
       hi_etaL2a->Fill(fabs(eta),anaW);	
       hi_ptres_L2->Fill((qL2/ptvL2-q/pt)/(q/pt));
       hi_chi2_L2->Fill(chi2L2);
       hi_prob_L2->Fill(probL2);
       if ( region == 0 ) hi_ptres_L2b->Fill((qL2/ptvL2-q/pt)/(q/pt));
       if ( region == 1 ) hi_ptres_L2o->Fill((qL2/ptvL2-q/pt)/(q/pt));
       if ( region == 2 ) hi_ptres_L2e->Fill((qL2/ptvL2-q/pt)/(q/pt));
       hi_ptpull_L2->Fill((1/ptvL2-1/pt)/EptL2); 
       //hi_thetapull_L2->Fill((t1.ThetavL2[idx_L2]-t1.Thetagen[idx_mu])/t1.EThetavL2[idx_L2]);
       //hi_phipull_L2->Fill((t1.PhivL2[idx_L2]-t1.Phigen[idx_mu])/t1.EPhivL2[idx_L2]);
     }
     
     if ( L3valid ) {
       hi_ptL3->Fill(pt,anaW);
       hi_etaL3->Fill(fabs(eta),anaW);
       hi_effL3->Fill(pt,eta,anaW);
       hi_ptL3a->Fill(pt,anaW);
       hi_etaL3a->Fill(fabs(eta),anaW);	
       hi_chi2_L3->Fill(chi2L3);
       hi_prob_L3->Fill(probL3);
       hi_ptres_L3->Fill((qL3/ptvL3-q/pt)/(q/pt));
       if ( region == 0 ) hi_ptres_L3b->Fill((qL3/ptvL3-q/pt)/(q/pt));
       if ( region == 1 ) hi_ptres_L3o->Fill((qL3/ptvL3-q/pt)/(q/pt));
       if ( region == 2 ) hi_ptres_L3e->Fill((qL3/ptvL3-q/pt)/(q/pt));
       hi_kink->Fill(kink,anaW/LumiLHC);
       hi_ptpull_L3->Fill((1/ptvL3-1/pt)/EptL3);
       //hi_thetapull_L3->Fill((t1.ThetavL3[idx_L3]-t1.Thetagen[idx_mu])/t1.EThetavL3[idx_L3]);
       //hi_phipull_L3->Fill((t1.PhivL3[idx_L3]-t1.Phigen[idx_mu])/t1.EPhivL3[idx_L3]);
     }
     
     // rates
     float step = 100./nbins;
     for ( int ibin = 0; ibin < nbins; ibin++ ) {
       float ptthres = ibin*step;
       if ( pt >= ptthres ) hi_ra->Fill(ptthres,anaW);
       if ( L2valid && ptL2 >= ptthres && ptL1 >= ptthres ) hi_l2->Fill(ptthres,anaW);
       if ( L3valid && ptL3 >= ptthres && ptL2 >= ptthres && ptL1 >= ptthres ) hi_l3->Fill(ptthres,anaW);
       
       
     }
     int nmax = hi_l1->GetNbinsX();
     for ( int ibin = 1; ibin <= nmax; ibin++ ) {
       float ptthres = hi_l1->GetBinLowEdge(ibin);
       if ( L1valid && ptL1 >= ptthres ) hi_l1->Fill(ptthres,anaW);
       if ( L1valid && ptL1 >= ptthres ) hi_l1a->Fill(ptthres,anaW);
     }  
     
   }
   
}



// ------------ method called once each job just before starting event loop  ------------
void 
RecoMuonEfficiency::beginJob(const edm::EventSetup&)
{
  cout << "DEBUG  240" << endl;
  online = false;cout << "DEBUG  241" << endl;

  effStyle = new TStyle("effStyle","Efficiency Study Style");
  
  //gROOT->Reset();cout << "DEBUG  243" << endl;
  effStyle->SetCanvasBorderMode(0);cout << "DEBUG  244" << endl;
  effStyle->SetCanvasBorderMode(0);cout << "DEBUG  245" << endl;
  effStyle->SetOptTitle(0);cout << "DEBUG  246" << endl;
  effStyle->SetOptStat(1000010);cout << "DEBUG  247" << endl;




  // Level-1 pt-scale (low edge)
  float ptvalues[32] = {
                           0.0,   1.5,   2.0,   2.5,   3.0,   3.5,
                           4.0,   4.5,   5.0,   6.0,   7.0,   8.0,
                          10.0,  12.0,  14.0,  16.0,  18.0,  20.0,
                          25.0,  30.0,  35.0,  40.0,  45.0,  50.0,
                          60.0,  70.0,  80.0,  90.0, 100.0, 120.0,
                         140.0, 200.0 };

 cout << "DEBUG  261" << endl;
  // efficiency 
  hi_ptgen  = new TH1F("hi_ptgen","gen pt",nbins,0.0,100.);
  hi_ptgen->Sumw2();
  hi_ptL1   = new TH1F("hi_ptL1","L1 pt",nbins,0.0,100.);
  hi_ptL1->Sumw2();
  hi_ptL2   = new TH1F("hi_ptL2","L2 pt",nbins,0.0,100.);
  hi_ptL2->Sumw2();
  hi_ptL3   = new TH1F("hi_ptL3","L3 pt",nbins,0.0,100.);
  hi_ptL3->Sumw2();
  hi_ptL2a  = new TH1F("hi_ptL2a","L2 pt",nbins,0.0,100.);
  hi_ptL2a->Sumw2();
  hi_ptL3a  = new TH1F("hi_ptL3a","L3 pt",nbins,0.0,100.);
  hi_ptL3a->Sumw2();
  cout << "DEBUG  275" << endl;
  hi_etagen = new TH1F("hi_etagen","gen eta",nbins,0,2.5);
  hi_etagen->Sumw2();
  hi_etaL1  = new TH1F("hi_etaL1","L1 eta",nbins,0,2.5);
  hi_etaL1->Sumw2();  
  hi_etaL2  = new TH1F("hi_etaL2","L2 eta",nbins,0,2.5);  
  hi_etaL2->Sumw2();
  hi_etaL3  = new TH1F("hi_etaL3","L3 eta",nbins,0,2.5);  
  hi_etaL3->Sumw2();
  hi_etaL2a = new TH1F("hi_etaL2a","L2 eta",nbins,0,2.5);
  hi_etaL2a->Sumw2();
  hi_etaL3a = new TH1F("hi_etaL3a","L3 eta",nbins,0,2.5);
  hi_etaL3a->Sumw2();
  cout << "DEBUG  288" << endl;
  hi_effgen = new TH2F("hi_effgen","gen",nbins,0.0,100.,nbins,-2.5,2.5);
  hi_effL1  = new TH2F("hi_effL1","L1",nbins,0.0,100.,nbins,-2.5,2.5);
  hi_effL2  = new TH2F("hi_effL2","L2",nbins,0.0,100.,nbins,-2.5,2.5);
  hi_effL3  = new TH2F("hi_effL3","L3",nbins,0.0,100.,nbins,-2.5,2.5);
  cout << "DEBUG  293" << endl;
  hi_ra = new TH1F("hi_ra","generator rates",nbins,0,100);
  hi_ra->Sumw2();  
  hi_l1 = new TH1F("hi_l1","L1 rates",31,ptvalues);
  hi_l1->Sumw2();
  hi_l1a = new TH1F("hi_l1a","L1 rates",nbins,0,100);
  hi_l1a->Sumw2();
  hi_l2 = new TH1F("hi_l2","L2 rates",nbins,0,100);
  hi_l2->Sumw2();
  hi_l3 = new TH1F("hi_l3","L3 rates",nbins,0,100);
  hi_l3->Sumw2();
  cout << "DEBUG  304" << endl;
  hi_l110pt = new TH1F("hi_l110pt","L1 10",nbins,0,100);
  hi_l110pt->Sumw2();
  hi_l120pt = new TH1F("hi_l120pt","L1 20",nbins,0,100);
  hi_l120pt->Sumw2();
  hi_l130pt = new TH1F("hi_l130pt","L1 30",nbins,0,100);
  hi_l130pt->Sumw2();
  hi_l140pt = new TH1F("hi_l140pt","L1 40",nbins,0,100);
  hi_l140pt->Sumw2();
  cout << "DEBUG  313" << endl;
  hi_l210pt = new TH1F("hi_l210pt","L2 10",nbins,0,100);
  hi_l210pt->Sumw2();
  hi_l220pt = new TH1F("hi_l220pt","L2 20",nbins,0,100);
  hi_l220pt->Sumw2();
  hi_l230pt = new TH1F("hi_l230pt","L2 30",nbins,0,100);
  hi_l230pt->Sumw2();
  hi_l240pt = new TH1F("hi_l240pt","L2 40",nbins,0,100);
  hi_l240pt->Sumw2();
  cout << "DEBUG  322" << endl;
  hi_l310pt = new TH1F("hi_l310pt","L3 10",nbins,0,100);
  hi_l310pt->Sumw2();
  hi_l320pt = new TH1F("hi_l320pt","L3 20",nbins,0,100);
  hi_l320pt->Sumw2();
  hi_l330pt = new TH1F("hi_l330pt","L3 30",nbins,0,100);
  hi_l330pt->Sumw2();
  hi_l340pt = new TH1F("hi_l340pt","L3 40",nbins,0,100);
  hi_l340pt->Sumw2();    
  cout << "DEBUG  331" << endl;
  // pt resolution
  hi_ptres_L1  = new TH1F("hi_ptres_L1","L1 pt resolution",nbins,-1.,1.);
  hi_ptres_L2  = new TH1F("hi_ptres_L2","L2 pt resolution",nbins,-1.,1.);
  hi_ptres_L3  = new TH1F("hi_ptres_L3","L3 pt resolution",nbins,-0.1,0.1);
  cout << "DEBUG  336" << endl;
  hi_ptres_L1b = new TH1F("hi_ptres_L1b","L1 pt resolution",nbins,-1.,1.);
  hi_ptres_L2b = new TH1F("hi_ptres_L2b","L2 pt resolution",nbins,-1.,1.);
  hi_ptres_L3b = new TH1F("hi_ptres_L3b","L3 pt resolution",nbins,-0.1,0.1);  
  hi_ptres_L1o = new TH1F("hi_ptres_L1o","L1 pt resolution",nbins,-1.,1.);
  hi_ptres_L2o = new TH1F("hi_ptres_L2o","L2 pt resolution",nbins,-1.,1.);
  hi_ptres_L3o = new TH1F("hi_ptres_L3o","L3 pt resolution",nbins,-0.1,0.1);   
  hi_ptres_L1e = new TH1F("hi_ptres_L1e","L1 pt resolution",nbins,-1.,1.);
  hi_ptres_L2e = new TH1F("hi_ptres_L2e","L2 pt resolution",nbins,-1.,1.);
  hi_ptres_L3e = new TH1F("hi_ptres_L3e","L3 pt resolution",nbins,-0.1,0.1);
  cout << "DEBUG  346" << endl;
  // chi2
  hi_chi2_L2 = new TH1F("hi_chi2_L2","Chi2 L2",nbins,0.,50.);
  hi_chi2_L3 = new TH1F("hi_chi2_L3","Chi2 L3",nbins,0.,20.);
  hi_prob_L2 = new TH1F("hi_prob_L2","Chi2 Prob L2",nbins,0.,0.0001);
  hi_prob_L3 = new TH1F("hi_prob_L3","Chi2 Prob L3",nbins,0.,1.); 
  cout << "DEBUG  352" << endl;
  // L2 pulls
  hi_ptpull_L2 = new TH1F("hi_ptpull_L2","L2 pt pull",nbins,-10.,10.);
  hi_thetapull_L2 = new TH1F("hi_thetapull_L2","L2 theta pull",nbins,-10.,10.);
  hi_phipull_L2 = new TH1F("hi_phipull_L2","L2 phi pull",nbins,-10.,10.);  
  cout << "DEBUG  357" << endl;
  // L3 pulls
  hi_ptpull_L3 = new TH1F("hi_ptpull_L3","L3 pt pull",nbins,-5.,5.);
  hi_thetapull_L3 = new TH1F("hi_thetapull_L3","L3 theta pull",nbins,-5.,5.);
  hi_phipull_L3 = new TH1F("hi_phipull_L3","L3 phi pull",nbins,-5.,5.);
  cout << "DEBUG  362" << endl;
  // kink
  hi_kink   = new TH1F("hi_kink","Kink L3",nbins,0.,50.);
  hi_regionR = new TH1F("hi_regionR","Region Of Interest #{Delta}R",100,0,10);
  cout << "DEBUG  365" << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RecoMuonEfficiency::endJob() {
  hi_ptL1 = (TH1F*)divideErr(hi_ptL1, hi_ptgen)->Clone("hi_ptL1");
  hi_ptL2 = (TH1F*)divideErr(hi_ptL2, hi_ptgen)->Clone("hi_ptL2");
  hi_ptL3 = (TH1F*)divideErr(hi_ptL3, hi_ptgen)->Clone("hi_ptL3");
 
  hi_ptL3a = divideErr(hi_ptL3a,hi_ptL2a);

  hi_etaL1 = (TH1F*)divideErr(hi_etaL1, hi_etagen)->Clone("hi_etaL1");
  hi_etaL2 = (TH1F*)divideErr(hi_etaL2, hi_etagen)->Clone("hi_etaL2");
  hi_etaL3 = (TH1F*)divideErr(hi_etaL3, hi_etagen)->Clone("hi_etaL3");  

  hi_etaL3a = divideErr(hi_etaL3a,hi_etaL2a);
  
  hi_effL1->Divide(hi_effgen);
  hi_effL2->Divide(hi_effgen);
  hi_effL3->Divide(hi_effgen);

  char* l1string = "Level-1";
  char* l2string = "Level-2";
  char* l3string = "Level-3";
  if ( !online ) {
    l2string = "Stand-alone";
    l3string = "Global";
  }

  //...................  
  effStyle->SetCanvasBorderMode(0);
  effStyle->SetPadBorderMode(1);
  effStyle->SetOptTitle(0);
  effStyle->SetStatFont(42);
  effStyle->SetTitleFont(22);
  effStyle->SetCanvasColor(10);
  effStyle->SetPadColor(0);
  effStyle->SetLabelFont(42,"x");
  effStyle->SetLabelFont(42,"y");
  effStyle->SetHistFillStyle(1001);
  effStyle->SetHistFillColor(0);
  effStyle->SetOptStat(0);
  //effStyle->SetOptStat(00011110);
  effStyle->SetOptFit(0111);
  effStyle->SetStatH(0.05); 
  //.................... 

  gROOT->SetStyle("effStyle");

  // 
  // plot histograms
  //
  TCanvas* c1 = new TCanvas("c1","Efficiency pt",10,10,700,500);
  c1->SetFillColor(0);
  c1->SetGrid(1);
  c1->SetTicky();
  c1->SetRightMargin(0.03);
  c1->SetTopMargin(0.02);
  c1->cd(); 
  
  hi_ptL1->SetXTitle("p_{T}^{#mu}");
  hi_ptL1->SetYTitle("Efficiency");
  hi_ptL1->SetTitleOffset(1.1,"x");
  hi_ptL1->SetTitleOffset(1.15,"y");
  hi_ptL1->SetMaximum(1.02);
  hi_ptL1->SetMinimum(0.5);
  hi_ptL1->SetStats(false);

  hi_ptL1->SetLineWidth(1);
  hi_ptL2->SetLineWidth(1);
  hi_ptL3->SetLineWidth(1);
    
  hi_ptL1->SetLineColor(2);
  hi_ptL2->SetLineColor(4);
  hi_ptL3->SetLineColor(1);

  hi_ptL1->SetFillColor(5);
  hi_ptL2->SetFillColor(7);
  hi_ptL3->SetFillColor(3);

  hi_ptL1->SetLineStyle(1);
  hi_ptL2->SetLineStyle(2);
  hi_ptL3->SetLineStyle(3);

  hi_ptL1 ->DrawCopy("HE");
  hi_ptL2 ->DrawCopy("HEsame");
  hi_ptL3 ->DrawCopy("HEsame");
  hi_ptL1 ->DrawCopy("AxisSame");

  TLegend* legend1 = new TLegend(0.6,0.2,0.8,0.4);
  legend1->SetTextAlign(32);
  legend1->SetTextColor(1);
  legend1->SetTextSize(0.04);
  legend1->AddEntry("hi_ptL1",l1string,"lf");
  legend1->AddEntry("hi_ptL2",l2string,"lf");
  legend1->AddEntry("hi_ptL3",l3string,"lf");
  legend1 ->Draw();
  c1->Write();
  // 
  TCanvas* c1a = new TCanvas("c1a","Algo Efficiency pt",10,10,700,500);
  c1a->SetFillColor(0);
  c1a->SetGrid(1);
  c1a->SetTicky();
  c1a->SetRightMargin(0.03);
  c1a->SetTopMargin(0.02);
  c1a->cd(); 
  
  hi_ptL3a->SetXTitle("p_{T}^{#mu}");
  hi_ptL3a->SetYTitle("Efficiency");
  hi_ptL3a->SetTitleOffset(1.1,"x");
  hi_ptL3a->SetStats(false);

  hi_ptL3a->SetMaximum(1.01);
  hi_ptL3a->SetMinimum(0.8);
  hi_ptL3a->SetLineColor(2);
  hi_ptL3a ->DrawCopy("HE");
  c1a->Write();
  //
  TCanvas* c2 = new TCanvas("c2","Efficiency eta",10,10,750,500);
  c2->SetFillColor(0);
  c2->SetGrid(1);
  c2->SetTicky();
  c2->SetRightMargin(0.03);
  c2->SetTopMargin(0.02);
  c2->cd(); 
  
  hi_etaL1->SetXTitle("|#eta^{#mu}|");
  hi_etaL1->SetYTitle("Efficiency"); 
  hi_etaL1->SetMaximum(1.01);
  hi_etaL1->SetMinimum(0.80);
  hi_etaL1->SetTitleOffset(1.1,"x");
  hi_etaL1->SetStats(false);

  hi_etaL1->SetLineWidth(1);
  hi_etaL2->SetLineWidth(1);
  hi_etaL3->SetLineWidth(1);

  hi_etaL1->SetLineColor(2);
  hi_etaL2->SetLineColor(4);
  hi_etaL3->SetLineColor(1);
  
  hi_etaL1->SetFillColor(5);
  hi_etaL2->SetFillColor(7);
  hi_etaL3->SetFillColor(3);
  
  hi_etaL1->SetLineStyle(1);
  hi_etaL2->SetLineStyle(2);
  hi_etaL3->SetLineStyle(3);
 
  hi_etaL1 ->DrawCopy("H");
  hi_etaL2 ->DrawCopy("Hsame");
  hi_etaL3 ->DrawCopy("Hsame");
  hi_etaL1 ->DrawCopy("AxisSame");
 
  TLegend* legend2 = new TLegend(0.6,0.2,0.8,0.4);
  legend2->SetTextAlign(32);
  legend2->SetTextColor(1);
  legend2->SetTextSize(0.04);
  legend2->AddEntry("hi_etaL1",l1string,"lf");
  legend2->AddEntry("hi_etaL2",l2string,"lf");
  legend2->AddEntry("hi_etaL3",l3string,"lf");
  legend2 ->Draw();
  c2->Write();
  // 
  TCanvas* c2a = new TCanvas("c2a","Algo Efficiency eta",10,10,700,500);
  c2a->SetFillColor(0);
  c2a->SetGrid(1);
  c2a->SetTicky();
  c2a->SetRightMargin(0.03);
  c2a->SetTopMargin(0.02);
  c2a->cd(); 
  
  hi_etaL3a->SetTitleOffset(1.1,"x");
  hi_etaL3a->SetStats(false);
  hi_etaL3a->SetXTitle("|#eta^{#mu}|");
  hi_etaL3a->SetYTitle("Efficiency"); 
  hi_etaL3a->SetMaximum(1.01);
  hi_etaL3a->SetMinimum(0.8);
  hi_etaL3a->SetLineColor(2);
  hi_etaL3a ->DrawCopy("H");
  c2a->Write();
/*  
  hi_ptres_L1;// ->Draw();
  hi_ptres_L2;// ->Draw();
  hi_ptres_L3;// ->Draw();
  
  hi_chi2_L2;// ->Draw();
  hi_chi2_L3;// ->Draw();
  
  hi_kink;// ->Draw();
*/

  //
  // turn on curves
  //
  TCanvas* c3 = new TCanvas("c3","turn on curves",20,20,720,465);
  c3->SetFillColor(0);
  c3->SetLogy(0);
  c3->SetGrid(1);
  c3->SetRightMargin(0.03);
  c3->SetTopMargin(0.02);
  c3->cd();

  TPad* pad1 = new TPad("pad1","Pad 1",0.0,0.5,0.5,1.0,10);
  TPad* pad2 = new TPad("pad2","Pad 2",0.5,0.5,1.0,1.0,10);
  TPad* pad3 = new TPad("pad3","Pad 3",0.0,0.0,0.5,0.5,10);
  TPad* pad4 = new TPad("pad4","Pad 4",0.5,0.0,1.0,0.5,10);
  pad1->SetRightMargin(0.03);
  pad1->SetTopMargin(0.02);
  pad1->SetGrid(1);
  pad1->SetLogy(0);
  pad2->SetRightMargin(0.03);
  pad2->SetTopMargin(0.02);
  pad2->SetGrid(1);
  pad2->SetLogy(0);
  pad3->SetRightMargin(0.03);
  pad3->SetTopMargin(0.02);
  pad3->SetGrid(1);
  pad3->SetLogy(0);
  pad4->SetRightMargin(0.03);
  pad4->SetTopMargin(0.02);
  pad4->SetGrid(1);
  pad4->SetLogy(0);
  pad1 ->Draw();
  pad2 ->Draw();
  pad3 ->Draw();
  pad4 ->Draw();

  hi_l110pt = divideErr(hi_l110pt, hi_ptgen);
  hi_l210pt = divideErr(hi_l210pt, hi_ptgen);
  hi_l310pt = divideErr(hi_l310pt, hi_ptgen);

  hi_l120pt = divideErr(hi_l120pt, hi_ptgen);
  hi_l220pt = divideErr(hi_l220pt, hi_ptgen);
  hi_l320pt = divideErr(hi_l320pt, hi_ptgen);

  hi_l130pt = divideErr(hi_l130pt, hi_ptgen);
  hi_l230pt = divideErr(hi_l230pt, hi_ptgen);
  hi_l330pt = divideErr(hi_l330pt, hi_ptgen);

  hi_l140pt = divideErr(hi_l140pt, hi_ptgen);
  hi_l240pt = divideErr(hi_l240pt, hi_ptgen);
  hi_l340pt = divideErr(hi_l340pt, hi_ptgen);

  pad1->cd();
  hi_l110pt->SetYTitle("Efficiency");
  hi_l110pt->SetXTitle("p_{T}^{#mu} [GeV/c]");
  hi_l110pt->GetXaxis()->SetRange(0,80);
  hi_l110pt->SetTitleOffset(1.2,"x");
  hi_l110pt->SetMinimum(0.0);
  hi_l110pt->SetMaximum(1.0);
  hi_l110pt->SetLineColor(2);
  hi_l110pt->SetLineStyle(1);
  hi_l110pt ->DrawCopy("HE");
  hi_l210pt->SetLineColor(4);
  hi_l210pt->SetLineStyle(2);
  hi_l210pt ->DrawCopy("HEsame");
  hi_l310pt->SetLineColor(8);
  hi_l310pt->SetLineStyle(3);
  hi_l310pt->SetLineWidth(3);
  hi_l310pt ->DrawCopy("HEsame");
  pad2->cd();
  hi_l120pt->SetYTitle("Efficiency");
  hi_l120pt->SetXTitle("p_{T}^{#mu} [GeV/c]");
  hi_l120pt->GetXaxis()->SetRange(0,80);
  hi_l120pt->SetTitleOffset(1.2,"x");
  hi_l120pt->SetMinimum(0.0);
  hi_l120pt->SetMaximum(1.0);
  hi_l120pt->SetLineColor(2);
  hi_l120pt->SetLineStyle(1);
  hi_l120pt ->DrawCopy("HE");
  hi_l220pt->SetLineColor(4);
  hi_l220pt->SetLineStyle(2);
  hi_l220pt ->DrawCopy("HEsame");
  hi_l320pt->SetLineColor(8);
  hi_l320pt->SetLineStyle(3);
  hi_l320pt->SetLineWidth(3);
  hi_l320pt ->DrawCopy("HEsame");
  pad3->cd();
  hi_l130pt->SetYTitle("Efficiency");
  hi_l130pt->SetXTitle("p_{T}^{#mu} [GeV/c]");
  hi_l130pt->GetXaxis()->SetRange(0,80);
  hi_l130pt->SetTitleOffset(1.2,"x");
  hi_l130pt->SetMinimum(0.0);
  hi_l130pt->SetMaximum(1.0);
  hi_l130pt->SetLineColor(2);
  hi_l130pt->SetLineStyle(1);
  hi_l130pt ->DrawCopy("HE");
  hi_l230pt->SetLineColor(4);
  hi_l230pt->SetLineStyle(2);
  hi_l230pt ->DrawCopy("HEsame");
  hi_l330pt->SetLineColor(8);
  hi_l330pt->SetLineStyle(3);
  hi_l330pt->SetLineWidth(3);
  hi_l330pt ->DrawCopy("HEsame");
  pad4->cd();
  hi_l140pt->SetYTitle("Efficiency");
  hi_l140pt->SetXTitle("p_{T}^{#mu} [GeV/c]");
  hi_l140pt->GetXaxis()->SetRange(0,80);
  hi_l140pt->SetTitleOffset(1.2,"x");
  hi_l140pt->SetMinimum(0.0);
  hi_l140pt->SetMaximum(1.0);
  hi_l140pt->SetLineColor(2);
  hi_l140pt->SetLineStyle(1);
  hi_l140pt ->DrawCopy("HE");
  hi_l240pt->SetLineColor(4);
  hi_l240pt->SetLineStyle(2);
  hi_l240pt ->DrawCopy("HEsame");
  hi_l340pt->SetLineColor(8);
  hi_l340pt->SetLineStyle(3);
  hi_l340pt->SetLineWidth(3);
  hi_l340pt ->DrawCopy("HEsame"); 
  c3->Write();

  //
  // resolution
  //
  TCanvas* c4 = new TCanvas("c4","1/pt Resolution L2",10,10,800,400);
  c4->SetFillColor(0);
  c4->SetGrid(1);
  c4->SetTicky();
  c4->SetRightMargin(0.03);
  c4->SetTopMargin(0.02);
  c4->cd();

  TPad* pad4a = new TPad("pad4a","Pad 1",0.00,0.0,0.33,1.0,10);
  TPad* pad4b = new TPad("pad4b","Pad 2",0.33,0.0,0.66,1.0,10);
  TPad* pad4c = new TPad("pad4c","Pad 3",0.66,0.0,0.99,1.0,10);
  pad4a ->Draw();
  pad4b ->Draw();
  pad4c ->Draw();

  pad4a->cd();
  TF1* g1 = new TF1("g1","gaus",-0.3,0.3);
  g1->SetLineColor(2);
  hi_ptres_L2b->Fit("g1","R");
  pad4b->cd();
  g1 = new TF1("g1","gaus",-0.4,0.45);
  g1->SetLineColor(2);
  hi_ptres_L2o->Fit("g1","R");
  pad4c->cd();
  g1 = new TF1("g1","gaus",-0.7,0.56);
  g1->SetLineColor(2);
  hi_ptres_L2e->Fit("g1","R"); 
  c4->Write();
  TCanvas* c5 = new TCanvas("c5","1/pt Resolution L3",10,10,800,400);
  c5->SetFillColor(0);
  c5->SetGrid(1);
  c5->SetTicky();
  c5->SetRightMargin(0.03);
  c5->SetTopMargin(0.02);
  c5->cd();

  TPad* pad5a = new TPad("pad5a","Pad 1",0.00,0.0,0.33,1.0,10);
  TPad* pad5b = new TPad("pad5b","Pad 2",0.33,0.0,0.66,1.0,10);
  TPad* pad5c = new TPad("pad5c","Pad 3",0.66,0.0,0.99,1.0,10);
  pad5a ->Draw();
  pad5b ->Draw();
  pad5c ->Draw();

  pad5a->cd();
  g1 = new TF1("g1","gaus",-0.05,0.05);
  g1->SetLineColor(2);
  hi_ptres_L3b->Fit("g1","R");
  pad5b->cd();
  g1 = new TF1("g1","gaus",-0.05,0.05);
  g1->SetLineColor(2);
  hi_ptres_L3o->Fit("g1","R");
  pad5c->cd();
  g1 = new TF1("g1","gaus",-0.1,0.1);
  g1->SetLineColor(2);
  hi_ptres_L3e->Fit("g1","R"); 
  c5->Write();

  //
  // pulls
  //
  TCanvas* c6 = new TCanvas("c6","L2 pulls",10,10,800,400);
  c6->SetFillColor(0);
  c6->SetGrid(1);
  c6->SetTicky();
  c6->SetRightMargin(0.03);
  c6->SetTopMargin(0.02);
  c6->cd();

  TPad* pad6a = new TPad("pad6a","Pad 1",0.00,0.0,0.33,1.0,10);
  TPad* pad6b = new TPad("pad6b","Pad 2",0.33,0.0,0.66,1.0,10);
  TPad* pad6c = new TPad("pad6c","Pad 3",0.66,0.0,0.99,1.0,10);
  pad6a ->Draw();
  pad6b ->Draw();
  pad6c ->Draw();

  pad6a->cd();
  g1 = new TF1("g1","gaus",-7.0,7.0);
  g1->SetLineColor(2);
  hi_ptpull_L2->Fit("g1","R");
  pad6b->cd();
  g1 = new TF1("g1","gaus",-5.0,5.0);
  g1->SetLineColor(2);
  hi_thetapull_L2->Fit("g1","R");
  pad6c->cd();
  g1 = new TF1("g1","gaus",-5.0,5.0);
  g1->SetLineColor(2);
  hi_phipull_L2->Fit("g1","R"); 
  c6->Write();
  TCanvas* c7 = new TCanvas("c7","L3 pulls",10,10,800,400);
  c7->SetFillColor(0);
  c7->SetGrid(1);
  c7->SetTicky();
  c7->SetRightMargin(0.03);
  c7->SetTopMargin(0.02);
  c7->cd();

  TPad* pad7a = new TPad("pad7a","Pad 1",0.00,0.0,0.33,1.0,10);
  TPad* pad7b = new TPad("pad7b","Pad 2",0.33,0.0,0.66,1.0,10);
  TPad* pad7c = new TPad("pad7c","Pad 3",0.66,0.0,0.99,1.0,10);
  pad7a ->Draw();
  pad7b ->Draw();
  pad7c ->Draw();

  pad7a->cd();
  g1 = new TF1("g1","gaus",-7.0,7.0);
  g1->SetLineColor(2);
  hi_ptpull_L3->Fit("g1","R");
  pad7b->cd();
  g1 = new TF1("g1","gaus",-5.0,5.0);
  g1->SetLineColor(2);
  hi_thetapull_L3->Fit("g1","R");
  pad7c->cd();
  g1 = new TF1("g1","gaus",-5.0,5.0);
  g1->SetLineColor(2);
  hi_phipull_L3->Fit("g1","R");
  c7->Write();
/*
  //
  // L1, L2, L3 rates
  //
  TCanvas* c10 = new TCanvas("c10","rates",20,20,720,465);
  c10->SetFillColor(0);
  c10->SetLogy(1);
  c10->SetGrid(1);
  c10->SetTicky();
  c10->SetRightMargin(0.03);
  c10->SetTopMargin(0.04);
  c10->cd();

  hi_ra->SetXTitle("p_{T}^{#mu} threshold [GeV/c]");
  hi_ra->SetYTitle("Rate [Hz]");
  hi_ra->SetTitleOffset(1.1,"x");
  hi_ra->SetTitleOffset(1.15,"y");
//  hi_ra->GetXaxis()->SetRange(0,80.0*nbins/100);
  hi_ra->SetMaximum(5000000);
  hi_ra->SetMinimum(0.1);
  hi_ra->SetLineColor(1);
  hi_ra;// ->Draw("H");

  hi_l1->SetMarkerStyle(20);
  hi_l1->SetMarkerColor(2);
  hi_l1->SetMarkerSize(0.8);
  hi_l1;// ->Draw("PEsame");
 
  hi_l2->SetMarkerStyle(21);
  hi_l2->SetMarkerColor(4);
  hi_l2->SetMarkerSize(0.8);
  hi_l2;// ->Draw("PEsame");

  hi_l3->SetMarkerStyle(29);
  hi_l3->SetMarkerColor(8);
  hi_l3->SetMarkerSize(1.1);
  hi_l3;// ->Draw("PEsame");
  
  TLegend* legend3 = new TLegend(0.68,0.7,0.88,0.9);
  legend3->SetTextAlign(12);
  legend3->SetTextColor(1);
  legend3->SetTextSize(0.03);
  legend3->AddEntry(hi_ra,"generator","l");
  legend3->AddEntry(hi_l1,"L1","p");
  legend3->AddEntry(hi_l2,"L2","p");
  legend3->AddEntry(hi_l3,"L3","p");
  legend3;// ->Draw();   

*/

  // efficiency 
  delete  hi_ptgen  ;
  delete  hi_ptL1   ;
  delete  hi_ptL2   ;
  delete  hi_ptL3   ;
  delete  hi_ptL2a  ;
  delete  hi_ptL3a  ;
     
  delete  hi_etagen ;
  delete  hi_etaL1  ;
  delete  hi_etaL2  ;
  delete  hi_etaL3  ;
  delete  hi_etaL2a ;
  delete  hi_etaL3a ;

  delete  hi_effgen ;
  delete  hi_effL1  ;
  delete  hi_effL2  ;
  delete  hi_effL3  ;

  delete  hi_ra ;
  delete  hi_l1 ;
  delete  hi_l1a ;
  delete  hi_l2 ;
  delete  hi_l3 ;

  delete  hi_l110pt ;
  delete  hi_l120pt ;
  delete  hi_l130pt ;
  delete  hi_l140pt ;

  delete  hi_l210pt ;
  delete  hi_l220pt ;
  delete  hi_l230pt ;
  delete  hi_l240pt ;

  delete  hi_l310pt ;
  delete  hi_l320pt ;
  delete  hi_l330pt ;
  delete  hi_l340pt ;

  // pt resolution
  delete  hi_ptres_L1  ;
  delete  hi_ptres_L2  ;
  delete  hi_ptres_L3  ;
  
  delete  hi_ptres_L1b ;
  delete  hi_ptres_L2b ;
  delete  hi_ptres_L3b ;
  delete  hi_ptres_L1o ;
  delete  hi_ptres_L2o ;
  delete  hi_ptres_L3o ;
  delete  hi_ptres_L1e ;
  delete  hi_ptres_L2e ;
  delete  hi_ptres_L3e ;

  // chi2
  delete  hi_chi2_L2 ;
  delete  hi_chi2_L3 ;
  delete  hi_prob_L2 ;
  delete  hi_prob_L3 ;
  
  // L2 pulls
  delete  hi_ptpull_L2 ;
  delete  hi_thetapull_L2 ;
  delete  hi_phipull_L2 ;

  // L3 pulls
  delete  hi_ptpull_L3 ;
  delete  hi_thetapull_L3 ;
  delete  hi_phipull_L3 ;

  // kink
  delete  hi_kink   ;
  hi_regionR->Write();
  delete hi_regionR;

  hFile->Write();
  hFile->Close();
  
}

//
// return h1/h2 with recalculated errors
//
TH1F* RecoMuonEfficiency::divideErr(TH1F* h1, TH1F* h2) {

  TH1F* hout = new TH1F(*h1);
  hout->Reset();
  hout->SetName("DivideErr");
  hout->Divide(h1,h2,1.,1.,"B");

  for (int i = 0; i <= hout->GetNbinsX()+1; i++ ) {
    Float_t tot   = h2->GetBinContent(i) ;
    Float_t tot_e = h2->GetBinError(i);
    Float_t eff = hout->GetBinContent(i) ;
    Float_t Err = 0.;
    if (tot > 0) Err = tot_e / tot * sqrt( eff* (1-eff) );
    if (eff == 1.) Err=1.e-3;
    hout->SetBinError(i, Err);
  }
  return hout;
}

float RecoMuonEfficiency::ptL2_90(float pt, float sigma_pt) {

  float s =  3.9; // 3.4
  return pt + s*sigma_pt*pt*pt;

}


float RecoMuonEfficiency::ptL3_90(float pt, float sigma_pt) {

  float s = 2.3; // 2.2
  return pt + s*sigma_pt*pt*pt;

}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoMuonEfficiency)
