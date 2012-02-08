#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

using namespace std;
using namespace reco;

class PatTest: public edm::EDAnalyzer 
{
  public:
    explicit PatTest(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~PatTest();

  private:
    edm::InputTag srcJets_,srcRho_,srcMet_;
};

PatTest::PatTest(edm::ParameterSet const& cfg) 
{
  srcJets_ = cfg.getParameter<edm::InputTag> ("jets");
  srcRho_  = cfg.getParameter<edm::InputTag> ("rho");
  srcMet_  = cfg.getParameter<edm::InputTag> ("met"); 
}  

void PatTest::beginJob() 
{
}

void PatTest::endJob() 
{
}

void PatTest::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(srcJets_,jets);
  edm::View<pat::Jet> pat_jets = *jets;

  edm::Handle<edm::View<MET> >  met;
  iEvent.getByLabel(srcMet_,met);

  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("goodOfflinePrimaryVertices",recVtxs);

  cout<<"nVtx = "<<recVtxs->size()<<", rho = "<<*rho<<", metSig = "<<(*met)[0].et()/(*met)[0].sumEt()<<endl;

  // ------ jet loop -----------------------------
  for(edm::View<pat::Jet>::const_iterator ijet = pat_jets.begin();ijet != pat_jets.end(); ++ijet) {
    float beta   = ijet->userFloat("beta");
    float ptD    = ijet->userFloat("ptD");
    float ptMax  = ijet->userFloat("ptMax");
    float axis1  = ijet->userFloat("axis1");
    float axis2  = ijet->userFloat("axis2");
    float tana   = ijet->userFloat("tana");
    float ttheta = ijet->userFloat("ttheta");
    cout<<"pt = "<<ijet->pt()<<", eta = "<<ijet->eta()<<", beta = "<<beta<<", ptD = "<<ptD<<", ptMax = "<<ptMax<<", axis1 = "<<axis1<<", axis2 = "<<axis2<<", tana = "<<tana<<", ttheta = "<<ttheta<<endl;
  }
}

PatTest::~PatTest() 
{
  
}

DEFINE_FWK_MODULE(PatTest);

