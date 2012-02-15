// -*- C++ -*-
//
// Package:    PatTracksOutOfJets
// Class:      PatTracksOutOfJets
// 
/**\class PatTracksOutOfJets PatTracksOutOfJets.cc Hbb/PatTracksOutOfJets/src/PatTracksOutOfJets.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesco FIORI
//         Created:  Thu Dec  1 16:01:56 CET 2011
// $Id: PatTracksOutOfJets.cc,v 1.1 2011/12/05 15:09:34 fiori Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//

// system include files
#include <memory>
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/HLTPrescaleTable.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

// user include files

#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include <DataFormats/PatCandidates/interface/Jet.h> 
#include <DataFormats/TrackReco/interface/Track.h> 
#include <DataFormats/MuonReco/interface/Muon.h> 
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Math/Point3D.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedRefCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedRefCandidateFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <TMath.h>
#include <TVector3.h>
//#include <TH1.h>
#include <TROOT.h>
#include <TFile.h>
#include "TNamed.h"
#include <vector>
#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TMatrixD.h>
#include <TVectorD.h>

//
// class declaration
//

class PatTracksOutOfJets : public edm::EDProducer {
   public:
      explicit PatTracksOutOfJets(const edm::ParameterSet&);
      ~PatTracksOutOfJets();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      edm::InputTag srcJets_ , srcVtx_ , srcTrk_ ;
      std::string   srcBtag_ ;


      // ----------member data ---------------------------
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
PatTracksOutOfJets::PatTracksOutOfJets(const edm::ParameterSet& iConfig)
{

  //produces<reco::TrackCollection>("outTracks");
  produces<reco::RecoChargedRefCandidateCollection>();

  srcJets_ = iConfig.getUntrackedParameter<edm::InputTag>  ("jets",edm::InputTag("jetExtender","extendedPatJets")) ;
  srcVtx_  = iConfig.getUntrackedParameter<edm::InputTag>  ("vtx",edm::InputTag("goodOfflinePrimaryVertices")) ;
  srcTrk_  = iConfig.getUntrackedParameter<edm::InputTag>  ("tracks",edm::InputTag("generalTracks")) ; 
  srcBtag_ = iConfig.getUntrackedParameter<std::string>    ("btagger",std::string("combinedSecondaryVertexBJetTags"));   

   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


PatTracksOutOfJets::~PatTracksOutOfJets()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PatTracksOutOfJets::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   //TrackCollection outTracks;
   std::auto_ptr<RecoChargedRefCandidateCollection> outTracks(new RecoChargedRefCandidateCollection);


   edm::Handle<edm::View<pat::Jet> > jets;
   iEvent.getByLabel(srcJets_,jets);
   edm::View<pat::Jet> pat_jets = *jets;

/*
   Handle<reco::JetTagCollection> PFbTagHandle3;
   iEvent.getByLabel("standardJetProbabilityak5PFBJetTags", PFbTagHandle3);
   const reco::JetTagCollection & PFbTags3 = *(PFbTagHandle3.product());
*/

   Handle<reco::VertexCollection> privtxs;
   int PriVtx=0;
   iEvent.getByLabel(srcVtx_,privtxs);
   const reco::VertexCollection & Pvtx= *(privtxs.product());
   PriVtx=Pvtx.size();

/*  string JetAlgorithmPF="ak5PFJets";
    Handle<PFJetCollection> jetHandle;
    iEvent.getByLabel(JetAlgorithmPF, jetHandle);
    const reco::PFJetCollection & jets = *(jetHandle.product());
*/
  
  reco::TrackRefVector pftracksin;
  pftracksin.clear();
  
  vector<TLorentzVector> vJetsAll;
  vJetsAll.clear();
  int JetsPF=0;
  vector<double> PFBJBP;
  
  double m=0.139;

  TLorentzVector a(0.0,0.0,0.0,0.0);
    
  vJetsAll.reserve(pat_jets.size());
  PFBJBP.reserve(pat_jets.size());

  for (unsigned int i = 0; i != pat_jets.size(); ++i) {
    
    a.SetPxPyPzE(pat_jets[i].px(),pat_jets[i].py(),pat_jets[i].pz(),pat_jets[i].energy());
    vJetsAll.push_back(a);
    
    reco::TrackRefVector pftrack = pat_jets[i].associatedTracks();  
    
    for(unsigned int pftr= 0;pftr<pftrack.size(); ++pftr) {
      if (i<4) pftracksin.push_back(pftrack[pftr]);
    }
    
    PFBJBP.push_back( pat_jets[i].bDiscriminator(srcBtag_) );
    ++JetsPF;
  }

  unsigned int track_activity= 0;
  float track_activity_pt = 0;
  
  if (JetsPF>3) {
    
     int ibtag[4];

     TMath::Sort(4,&(PFBJBP.at(0)),ibtag);
    
     TLorentzVector *bjet0 = &vJetsAll[ibtag[0]];
     TLorentzVector *bjet1 = &vJetsAll[ibtag[1]];

     TLorentzVector t(0.0,0.0,0.0,0.0);
     Handle<TrackCollection> tracks;
     iEvent.getByLabel(srcTrk_, tracks);

     // get the primary vertices points
     std::vector<math::XYZPoint> vpoints;
     vpoints.reserve(Pvtx.size());
     for (VertexCollection::const_iterator vtx=Pvtx.begin();vtx!=Pvtx.end();++vtx) {
       vpoints.push_back(vtx->position());       
     }

     int j=0;
     for(TrackCollection::const_iterator itTrack = tracks->begin();itTrack != tracks->end(); ++itTrack,++j){
       bool out=true;
       // only high purity tracks with pT > 300 MeV
       if (itTrack->quality(TrackBase::highPurity) && itTrack->pt()>0.3){
	 for (unsigned int i1 = 0; i1 != pftracksin.size(); ++i1) {
	   if (&*itTrack== &*(pftracksin[i1]))  out=false;
	 }
	 
	 float dzPV0 = itTrack->dz(vpoints[0]);
	 float dzErr = itTrack->dzError();
	 // minimum z-distance of track to the first PV : 1mm && 3sigma
	 if (fabs(dzPV0) > 0.2 || fabs(dzPV0/dzErr)>3) out = false;
	 // loop over secondary (softer) Primary Vertices and exclude tracks more compatible with those 

	 for (int iv =1; iv <PriVtx; ++iv) {
	   float dz =  itTrack->dz(vpoints[iv]);
	   if (fabs(dz)<fabs(dzPV0)) out = false;
	 }

	 float dR0 = 0.5;
	 double energy=sqrt(0.139*0.139+itTrack->px()*itTrack->px()+itTrack->py()*itTrack->py()+itTrack->pz()*itTrack->pz());
	 t.SetPxPyPzE(itTrack->px(),itTrack->py(),itTrack->pz(),energy);
	 // remove tracks in an ellipse area around the two b-jets
	 float dRbb = bjet0->DeltaR(*bjet1);
	 float dRbjet0 = t.DeltaR(*bjet0); 
	 float dRbjet1 = t.DeltaR(*bjet1); 
	 if (dRbjet0 + dRbjet1 < dRbb + 2*dR0) out = false;
	 if (out) {
	   // tracks with pT > 1 GeV go in the track activity variables 
	   if (itTrack->pt() > 1) {
	     track_activity_pt += itTrack->pt();
	   }
	   track_activity++;
	   TrackRef Tref=TrackRef(tracks,j);
	   RecoChargedRefCandidate refCand=RecoChargedRefCandidate(Tref,m);
	   outTracks->push_back(refCand);
	 }
       }
     }
  }


  /*RecoChargedRefCandidateCollection outCand;
  double m=0.139;
  for (int k=0; k<outTracks->size(); k++){
    outCand.push_back(RecoChargedRefCandidate(outTracks[k].trackRef(),m));
    }*/
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(pOut);
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
  //std::auto_ptr<TrackCollection> pOut(new TrackCollection(*outTracks));
  //iEvent.put(outTracks);
  iEvent.put(outTracks);
}

// ------------ method called once each job just before starting event loop  ------------
void 
PatTracksOutOfJets::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PatTracksOutOfJets::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
PatTracksOutOfJets::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PatTracksOutOfJets::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PatTracksOutOfJets::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PatTracksOutOfJets::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PatTracksOutOfJets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PatTracksOutOfJets);
