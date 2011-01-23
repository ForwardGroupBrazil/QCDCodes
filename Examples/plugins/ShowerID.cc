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
#include "DataFormats/MuonReco/interface/MuonShower.h"


class ShowerID : public edm::EDProducer {
   public:
      explicit ShowerID(const edm::ParameterSet&);
      ~ShowerID();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
edm::InputTag src_;
std::string result_;
};


ShowerID::ShowerID(const edm::ParameterSet& iConfig)
{
src_= iConfig.getParameter<edm::InputTag>("src");
result_ = iConfig.getParameter<std::string>("result");
produces<edm::ValueMap<float> >().setBranchAlias("ShowerInformation");

}


ShowerID::~ShowerID()
{
 
}


// ------------ method called to produce the data  ------------
void
ShowerID::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  Handle<edm::ValueMap<reco::MuonShower> > ShowerMap;
  iEvent.getByLabel( src_, ShowerMap );
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel("muons",muons);
  std::vector<float> values;
  values.reserve(muons->size());
  
  unsigned int muonIdx = 0;
  for(reco::MuonCollection::const_iterator muon = muons->begin();
      muon != muons->end(); ++muon) {
    reco::MuonRef muonRef(muons, muonIdx);
    reco::MuonShower muonShowerInformation = (*ShowerMap)[muonRef];
    if(result_ == "allHits") {
         values.push_back((muonShowerInformation.nStationHits).at(0));
    }
    if(result_ == "correlatedHits") values.push_back((muonShowerInformation.nStationCorrelatedHits).at(0));
    if(result_ == "showerSize") values.push_back((muonShowerInformation.stationShowerSizeT).at(0));
    if(result_ == "showerDeltaR") values.push_back((muonShowerInformation.stationShowerDeltaR).at(0));
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
ShowerID::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ShowerID::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ShowerID);
