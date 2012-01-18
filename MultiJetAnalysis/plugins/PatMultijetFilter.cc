#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

using namespace std;
using namespace reco;
using namespace edm;

class PatMultijetFilter : public edm::EDFilter 
{
  public:
    typedef std::vector<Jet> JetCollection; 
    explicit PatMultijetFilter(const edm::ParameterSet & cfg);
    virtual void beginJob();
    bool filter(edm::Event &event, const edm::EventSetup &iSetup);
    virtual void endJob();  
  private:
    edm::InputTag jets_;
    double minPt_;
    int minNjets_;
};

PatMultijetFilter::PatMultijetFilter(ParameterSet const& cfg) 
{
  jets_     = cfg.getParameter<edm::InputTag>("jets");
  minNjets_ = cfg.getParameter<int>("minNjets");
  minPt_    = cfg.getParameter<double>("minPt");
}
////////////////////////////////////////////////////////////////////
void PatMultijetFilter::beginJob() 
{
}
////////////////////////////////////////////////////////////////////
bool PatMultijetFilter::filter(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jets_,jets);
  edm::View<pat::Jet> pat_jets = *jets;

  edm::View<pat::Jet>::const_iterator i_jet;
  int N(0);
  bool check(false);
  //////// Find the number of real jets (pt>ptTrhreshold) /////////// 
  for(i_jet = pat_jets.begin(); i_jet != pat_jets.end(); i_jet++) {
    double pt = i_jet->pt();
    if (pt > minPt_) {
      N++;
    }
  }
  if (N >= minNjets_) {
    check = true; 
  } 
  return check;
}
//////////////////////////////////////////////////////////////////////////////////////////
void PatMultijetFilter::endJob() 
{
}
DEFINE_FWK_MODULE(PatMultijetFilter);
