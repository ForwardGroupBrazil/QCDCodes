// system include files
#include <memory>
#include <iostream>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"


class JetBetaProducer : public edm::EDProducer {
   public:
      explicit JetBetaProducer(const edm::ParameterSet&);
      ~JetBetaProducer();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data --------------------------
        edm::InputTag src_;
        std::string result_;
};

JetBetaProducer::JetBetaProducer(const edm::ParameterSet& iConfig)
{
        src_    = iConfig.getParameter<edm::InputTag>("jets");
        produces<edm::ValueMap<float> >().setBranchAlias("beta");
}

JetBetaProducer::~JetBetaProducer()
{

}

// ------------ method called to produce the data  ------------
void
JetBetaProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  
  edm::Handle<reco::PFJetCollection> pfjets;
  iEvent.getByLabel(src_,pfjets);
  std::vector<float> values;
  values.reserve(pfjets->size());
  //-------------- Vertex Info -----------------------------------
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("goodOfflinePrimaryVertices",recVtxs);
    
  unsigned int jetIdx = 0;
  for(reco::PFJetCollection::const_iterator ijet = pfjets->begin(); ijet != pfjets->end(); ++ijet) {
    //---- vertex association -----------
    //---- get the vector of tracks -----
    reco::TrackRefVector vTrks(ijet->getTrackRefs());
    float sumTrkPt(0.0),sumTrkPtBeta(0.0),beta(0.0);
    //---- loop over the tracks of the jet ----
    for(reco::TrackRefVector::const_iterator i_trk = vTrks.begin(); i_trk != vTrks.end(); i_trk++) {
      sumTrkPt += (*i_trk)->pt();
      //---- loop over all vertices ----------------------------
      for(unsigned ivtx = 0;ivtx < recVtxs->size();ivtx++) {
        //---- loop over the tracks associated with the vertex ---
        for(reco::Vertex::trackRef_iterator i_vtxTrk = (*recVtxs)[ivtx].tracks_begin(); i_vtxTrk != (*recVtxs)[ivtx].tracks_end(); ++i_vtxTrk) {
          //---- match the jet track to the track from the vertex ----
          reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
          //---- check if the tracks match -------------------------
          if (trkRef == (*i_trk)) {
            if (ivtx > 0) {
              sumTrkPtBeta += (*i_trk)->pt();
            }   
            break;
          }
        }
      } 
    }
    if (sumTrkPt > 0) {
      beta = 1.-sumTrkPtBeta/sumTrkPt;  
    }
    values.push_back(beta);
     ++jetIdx;
  }

  std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*out);
  filler.insert(pfjets, values.begin(), values.end());
  filler.fill();

  // put value map into event
  iEvent.put(out);
}

// ------------ method called once each job just before starting event loop  ------------
void
JetBetaProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
JetBetaProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetBetaProducer);
