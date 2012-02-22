#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

using namespace std;
using namespace edm;
using namespace reco;

class PartonFilter : public edm::EDFilter 
{
  public: 
    explicit PartonFilter(const edm::ParameterSet & cfg);
    virtual void beginJob();
    bool filter(edm::Event &event, const edm::EventSetup &iSetup);
    virtual void endJob();  
  private:
    edm::InputTag src_;
    int nPartons_,partonType_;
};

PartonFilter::PartonFilter(ParameterSet const& cfg) 
{
  src_        = cfg.getParameter<edm::InputTag>("src");
  nPartons_   = cfg.getParameter<int>          ("nPartons");
  partonType_ = cfg.getParameter<int>          ("partonType");
}
////////////////////////////////////////////////////////////////////
void PartonFilter::beginJob() 
{
}
////////////////////////////////////////////////////////////////////
bool PartonFilter::filter(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  edm::Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel(src_,genParticles);

  int N(0);
  bool check(false);
  for(unsigned ip = 0; ip < genParticles->size(); ++ ip) {
    const GenParticle &p = (*genParticles)[ip];
    if (p.status() != 3) continue;
    if (p.pdgId() != partonType_) continue;
    int ndst3(0);
    for(unsigned k = 0; k < p.numberOfDaughters(); k++) {
      if (p.daughter(k)->status() == 3) {
        ndst3++; 
      }
    }
    if (ndst3 > 0) continue;
    N++; 
  }
  if (N >= nPartons_) {
    check = true; 
  } 
  return check;
}
//////////////////////////////////////////////////////////////////////////////////////////
void PartonFilter::endJob() 
{
}
DEFINE_FWK_MODULE(PartonFilter);
