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

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"


class MuonCosmicMaps : public edm::EDProducer {
   public:
      explicit MuonCosmicMaps(const edm::ParameterSet&);
      ~MuonCosmicMaps();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
edm::InputTag src_;
std::string result_;
};


MuonCosmicMaps::MuonCosmicMaps(const edm::ParameterSet& iConfig)
{
src_= iConfig.getParameter<edm::InputTag>("src");
result_ = iConfig.getParameter<std::string>("result");
produces<edm::ValueMap<float> >().setBranchAlias("CosmicDiscriminators");

}


MuonCosmicMaps::~MuonCosmicMaps()
{
 
}


// ------------ method called to produce the data  ------------
void
MuonCosmicMaps::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  Handle<edm::ValueMap<reco::MuonCosmicCompatibility> > CosmicMap;
  iEvent.getByLabel( src_, CosmicMap );
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel("muons",muons);
  std::vector<float> values;
  values.reserve(muons->size());
  
  unsigned int muonIdx = 0;
  for(reco::MuonCollection::const_iterator muon = muons->begin();
      muon != muons->end(); ++muon) {
    reco::MuonRef muonRef(muons, muonIdx);
    reco::MuonCosmicCompatibility muonCosmicCompatibility = (*CosmicMap)[muonRef];
    if(result_ == "cosmicCompatibility") values.push_back(muonCosmicCompatibility.cosmicCompatibility);
    if(result_ == "timeCompatibility") values.push_back(muonCosmicCompatibility.timeCompatibility);
    if(result_ == "backToBackCompatibility") values.push_back(muonCosmicCompatibility.backToBackCompatibility);
    if(result_ == "overlapCompatibility") values.push_back(muonCosmicCompatibility.overlapCompatibility);
    //if(result_ == "ipCompatibility") values.push_back(muonCosmicCompatibility.ipCompatibility);
    //if(result_ == "vertexCompatibility") values.push_back(muonCosmicCompatibility.vertexCompatibility);
    ++muonIdx;
  }
  
  std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*out);
  filler.insert(muons, values.begin(), values.end());
  filler.fill();
  
  // put value map into event
  iEvent.put(out);
  
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonCosmicMaps::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonCosmicMaps::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonCosmicMaps);
