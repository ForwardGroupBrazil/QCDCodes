// -*- C++ -*-
//
// Package:    ParamSwitch
// Class:      ParamSwitch
// 
/**\class ParamSwitch ParamSwitch.cc UserCode/ParamSwitch/src/ParamSwitch.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam A Everett
//         Created:  Fri Dec 12 09:18:48 EST 2008
// $Id: ParamSwitch.cc,v 1.1 2009/01/05 21:39:35 aeverett Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

//
// class decleration
//

class ParamSwitch : public edm::EDAnalyzer {
   public:
      explicit ParamSwitch(const edm::ParameterSet&);
      ~ParamSwitch();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  edm::InputTag muonsTag;
  edm::InputTag selectionTag;
  const edm::ParameterSet parset_;   

  double chi2Cut_;
  double maxPtCut_;
  double minPtCut_;
  int hitCut_;

  bool hitFlag_;
  bool chiFlag_;
  bool ptFlag_;

  TH1F *h_switch_tracker;
  TH1F *h_switch_global;
  TH1F *h_switch_tracker_barrel;
  TH1F *h_switch_global_barrel;
  TH1F *h_switch_tracker_endcap;
  TH1F *h_switch_global_endcap;
  TH1F *h_switch_tracker_overlap;
  TH1F *h_switch_global_overlap;
  TH1F *h_switch_50;
  TH1F *h_switch_100;
  TH1F *h_switch_150;
  TH1F *h_switch_200;
  TH1F *h_switch_250;
  TH1F *h_switch_300;
  TH1F *h_switch2_50;
  TH1F *h_switch2_100;
  TH1F *h_switch2_150;
  TH1F *h_switch2_200;
  TH1F *h_switch2_250;
  TH1F *h_switch2_300;

  TH1F *h_switch_s1;
  TH1F *h_switch_s2;
  TH1F *h_switch_s3;

  TH1I *h_switch_q_s1;
  TH1I *h_switch_q_s2;
  TH1I *h_switch_q_s3;

  TH1F *h_switch_glbParam;

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
ParamSwitch::ParamSwitch(const edm::ParameterSet& iConfig) :
  muonsTag(    iConfig.getParameter< edm::InputTag >("muonsTag")),
  selectionTag(iConfig.getParameter< edm::InputTag >("selectionTag")),
  parset_(iConfig)

{
   //now do what ever initialization is needed
  //std::cout <<"l 100" << std::endl;
  maxPtCut_ =  iConfig.getParameter<double>("maxPtCut");
  minPtCut_ =  iConfig.getParameter<double>("minPtCut");
  hitCut_ =  iConfig.getParameter<int>("hitCut");
  chi2Cut_ =  iConfig.getParameter<double>("chi2Cut");
  chiFlag_ = iConfig.getParameter<bool>("chiFlag");
  hitFlag_ = iConfig.getParameter<bool>("hitFlag");
  ptFlag_ = iConfig.getParameter<bool>("ptFlag");
}


ParamSwitch::~ParamSwitch()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ParamSwitch::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  //cout <<"l 123" << endl;
  edm::Handle<reco::MuonCollection> muonCollectionH;
  
  edm::LogVerbatim("MuonTrackSelector") << std::endl << "muonsTag = " << muonsTag.label() << std::endl;
  edm::LogVerbatim("MuonTrackSelector") << "selectionTag = " << selectionTag.label() << std::endl;
  //cout <<"l 129" << endl;
  iEvent.getByLabel(muonsTag,muonCollectionH);
  if(muonCollectionH->size()==0) {
    edm::LogVerbatim("MuonTrackSelector") << "\nNo muons in this event\n";
  }
  //cout <<"l 134" << endl;
  for(reco::MuonCollection::const_iterator muon = muonCollectionH->begin(); muon != muonCollectionH->end(); ++muon) {
    bool isGoodSelection = false;
    bool passSpecial = true;
    //cout <<"l 137" << endl;
    if (selectionTag.label() == "AllGlobalMuons") {
      isGoodSelection = muon->isGood(reco::Muon::AllGlobalMuons);
      //cout <<"l 140" << endl;
      if(isGoodSelection) {
	//cout <<"l 142" << endl;
	TrackRef trkRef_inner = muon->innerTrack();//cout <<"l 143" << endl;
	TrackRef trkRef_outer = muon->outerTrack();//cout <<"l 144" << endl;
	TrackRef trkRef_global = muon->globalTrack();//cout <<"l 145" << endl;

	//////////////////
	//////////////////
	bool returnChiVal = false;
	double chiVal = trkRef_global->normalizedChi2();
	if (chi2Cut_ > 0.) { 
	  if(chiVal <= fabs(chi2Cut_) ) returnChiVal = true;
	} else {
	  if(chiVal >= fabs(chi2Cut_) ) returnChiVal = true;
	}
	//////////////////
	bool returnHitVal = false;
	int hitTk = muon->track().get()->recHitsSize();
	int hitGlbTk =   muon->combinedMuon().get()->hitPattern().numberOfValidTrackerHits();
	int hitSta =  muon->standAloneMuon().get()->recHitsSize();
	int hitGlbSta = muon->combinedMuon().get()->hitPattern().numberOfValidMuonHits();
	int hitGlb =  muon->combinedMuon().get()->hitPattern().numberOfValidHits();
	
	int missingSta = hitSta-hitGlbSta;
	double fractionSta = (hitSta > 0) ? (hitSta-hitGlbSta)/hitSta : 0;
	int missingTk = hitTk-hitGlbTk;
	double fractionTk = (hitTk > 0) ? (hitTk-hitGlbTk)/hitTk : 0;
	
	if(hitCut_ > 0) {
	  if(missingSta >= hitCut_ || missingTk >= hitCut_) returnHitVal = true;
	} else {
	  if(missingSta >= hitCut_ || missingTk >= hitCut_) returnHitVal = true;
	}
	//////////////////
	bool returnPtVal = false;
	if(muon->pt() > minPtCut_ && muon->pt() < maxPtCut_) returnPtVal = true;
	//////////////////
	//////////////////

	passSpecial = ( (hitFlag_ && returnHitVal) ||
			(chiFlag_ && returnChiVal) ||
			(ptFlag_  && returnPtVal ) ) ? true : false;

	passSpecial = ( !hitFlag_ && !chiFlag_ && !ptFlag_) ? true : false;

	if (passSpecial == false) continue;

	//cout <<"l 146" << endl;
	float innerPt = trkRef_inner->pt();
	float globalPt = trkRef_global->pt();
	float globalEta =  trkRef_global->eta();
	//cout <<"l 149" << endl;
	h_switch_tracker->Fill(innerPt);
	h_switch_global->Fill(globalPt);
	if(globalEta <= 0.8) {
	  h_switch_tracker_barrel->Fill(innerPt);
	  h_switch_global_barrel->Fill(globalPt);
	} else if (globalEta > 0.8 && globalEta <= 1.2) {
	  h_switch_tracker_overlap->Fill(innerPt);
	  h_switch_global_overlap->Fill(globalPt);
	} else if (globalEta > 1.2 && globalEta <= 2.4) {
	  h_switch_tracker_endcap->Fill(innerPt);
	  h_switch_global_endcap->Fill(globalPt);
	}
	//cout <<"l 152" << endl;
	float fillPt50  = innerPt < 50  ? innerPt : globalPt;
	float fillPt100 = innerPt < 100 ? innerPt : globalPt;
	float fillPt150 = innerPt < 150 ? innerPt : globalPt;
	float fillPt200 = innerPt < 200 ? innerPt : globalPt;
	float fillPt250 = innerPt < 250 ? innerPt : globalPt;
	float fillPt300 = innerPt < 300 ? innerPt : globalPt;
	float fillPt50_2  = innerPt >  50 && globalPt >  50 ? globalPt : innerPt;
	float fillPt100_2 = innerPt > 100 && globalPt > 100 ? globalPt : innerPt;
	float fillPt150_2 = innerPt > 150 && globalPt > 150 ? globalPt : innerPt;
	float fillPt200_2 = innerPt > 200 && globalPt > 200 ? globalPt : innerPt;
	float fillPt250_2 = innerPt > 250 && globalPt > 250 ? globalPt : innerPt;
	float fillPt300_2 = innerPt > 300 && globalPt > 300 ? globalPt : innerPt;
	//cout <<"l 159" << endl;
	h_switch_50->Fill(fillPt50);
	h_switch_100->Fill(fillPt100);
	h_switch_150->Fill(fillPt150);
	h_switch_200->Fill(fillPt200);
	h_switch_250->Fill(fillPt250);
	h_switch_300->Fill(fillPt300);
	h_switch2_50->Fill(fillPt50_2);
	h_switch2_100->Fill(fillPt100_2);
	h_switch2_150->Fill(fillPt150_2);
	h_switch2_200->Fill(fillPt200_2);
	h_switch2_250->Fill(fillPt250_2);
	h_switch2_300->Fill(fillPt300_2);
	//cout <<"l 166" << endl;
	float sigma = trkRef_inner->qoverpError();
	float innerQoverP = trkRef_inner->qoverp(); 
	float globalQoverP = trkRef_global->qoverp();
	float fillPtS1 = (fabs( innerQoverP-globalQoverP ) > sigma) ? innerPt : globalPt;
	float fillPtS2 = (fabs( innerQoverP-globalQoverP ) > 2*sigma) ? innerPt : globalPt;
	float fillPtS3 = (fabs( innerQoverP-globalQoverP ) > 3*sigma) ? innerPt : globalPt;
	h_switch_s1->Fill(fillPtS1);
	h_switch_s2->Fill(fillPtS2);
	h_switch_s3->Fill(fillPtS3);

	int global_q = trkRef_global->charge();
	int inner_q = trkRef_inner->charge();
	int fillQS1 = (fabs( innerQoverP-globalQoverP ) > sigma) ? inner_q : global_q;
	int fillQS2 = (fabs( innerQoverP-globalQoverP ) > 2*sigma) ? inner_q : global_q;
	int fillQS3 = (fabs( innerQoverP-globalQoverP ) > 3*sigma) ? inner_q : global_q;
	h_switch_q_s1->Fill(fillQS1*inner_q);
	h_switch_q_s2->Fill(fillQS2*inner_q);
	h_switch_q_s3->Fill(fillQS3*inner_q);

	float fillGlbParam = (missingSta >= 5 || missingTk >= 5 || chiVal >= 2.5) ? innerPt : globalPt ;
	h_switch_glbParam->Fill(fillGlbParam);
      }
    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
ParamSwitch::beginJob(const edm::EventSetup&)
{
 edm::Service<TFileService> fs;

 h_switch_tracker = fs->make<TH1F>("h_switch_tracker","Pt spectrum for tracker",700,0,700);
 h_switch_global = fs->make<TH1F>("h_switch_global","Pt spectrum for global",700,0,700);

 h_switch_tracker_barrel = fs->make<TH1F>("h_switch_tracker_barrel","Pt spectrum for tracker (barrel)",700,0,700);
 h_switch_global_barrel = fs->make<TH1F>("h_switch_global_barrel","Pt spectrum for global (barrel)",700,0,700);

 h_switch_tracker_overlap = fs->make<TH1F>("h_switch_tracker_overlap","Pt spectrum for tracker (overlap)",700,0,700);
 h_switch_global_overlap = fs->make<TH1F>("h_switch_global_overlap","Pt spectrum for global (overlap)",700,0,700);

 h_switch_tracker_endcap = fs->make<TH1F>("h_switch_tracker_endcap","Pt spectrum for tracker (endcap)",700,0,700);
 h_switch_global_endcap = fs->make<TH1F>("h_switch_global_endcap","Pt spectrum for global (endcap)",700,0,700);

 h_switch_50 = fs->make<TH1F>("h_switch_50","Pt spectrum with switch at 50 GeV",700,0,700);
 h_switch_100 = fs->make<TH1F>("h_switch_100","Pt spectrum with switch at 100 GeV",700,0,700);
 h_switch_150 = fs->make<TH1F>("h_switch_150","Pt spectrum with switch at 150 GeV",700,0,700);
 h_switch_200 = fs->make<TH1F>("h_switch_200","Pt spectrum with switch at 200 GeV",700,0,700);
 h_switch_250 = fs->make<TH1F>("h_switch_250","Pt spectrum with switch at 250 GeV",700,0,700);
 h_switch_300 = fs->make<TH1F>("h_switch_300","Pt spectrum with switch at 300 GeV",700,0,700);

 h_switch2_50 = fs->make<TH1F>("h_switch2_50","Pt spectrum with switch_{2} at 50 GeV",700,0,700);
 h_switch2_100 = fs->make<TH1F>("h_switch2_100","Pt spectrum with switch_{2} at 100 GeV",700,0,700);
 h_switch2_150 = fs->make<TH1F>("h_switch2_150","Pt spectrum with switch_{2} at 150 GeV",700,0,700);
 h_switch2_200 = fs->make<TH1F>("h_switch2_200","Pt spectrum with switch_{2} at 200 GeV",700,0,700);
 h_switch2_250 = fs->make<TH1F>("h_switch2_250","Pt spectrum with switch_{2} at 250 GeV",700,0,700);
 h_switch2_300 = fs->make<TH1F>("h_switch2_300","Pt spectrum with switch_{2} at 300 GeV",700,0,700);

 h_switch_s1 = fs->make<TH1F>("h_switch_s1","Pt spectrum with switch at 1 sigma",700,0,700);
 h_switch_s2 = fs->make<TH1F>("h_switch_s2","Pt spectrum with switch at 2 sigma",700,0,700);
 h_switch_s3 = fs->make<TH1F>("h_switch_s3","Pt spectrum with switch at 3 sigma",700,0,700);

 h_switch_glbParam = fs->make<TH1F>("h_switch_glbParam","Pt spectrum with global parameter switch",700,0,700);

 h_switch_q_s1 = fs->make<TH1I>("h_switch_q_s1","q_{default} * q_{Tk} with switch at 1 sigma",5,-2,2);
 h_switch_q_s2 = fs->make<TH1I>("h_switch_q_s2","q_{default} * q_{Tk} with switch at 2 sigma",5,-2,2);
 h_switch_q_s3 = fs->make<TH1I>("h_switch_q_s3","q_{default} * q_{Tk} with switch at 3 sigma",5,-2,2);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ParamSwitch::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ParamSwitch);
