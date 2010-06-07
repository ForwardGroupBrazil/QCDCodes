// -*- C++ -*-
//
// Package:    MuonHitCounter
// Class:      MuonHitCounter
// 
/**\class MuonHitCounter MuonHitCounter.cc PhysicsTools/PatAlgos/src/MuonHitCounter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Nov 16 16:12 (lxplus231.cern.ch)
//         Created:  Sun Nov 16 16:14:09 CET 2008
// $Id: MuonHitCounter.cc,v 1.1 2010/06/06 17:26:02 aeverett Exp $
//
//


// system include files
#include <memory>
#include <set>
#include <ext/hash_map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "UserCode/Examples/interface/muonHitCount.h"

#include "boost/lexical_cast.hpp"

using boost::lexical_cast;

//
// class decleration
class MuonHitCounter : public edm::EDProducer {
    public:
        explicit MuonHitCounter(const edm::ParameterSet&);
        ~MuonHitCounter();

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);

        /// Write a ValueMap<int> in the event
        void writeValueMap(edm::Event &iEvent,
                const edm::Handle<edm::View<reco::Muon> > & handle,
                const std::vector<int> & values,
                const std::string    & label) const ;

        /// The muons
        edm::InputTag src_;

        /// Use global track instead of standalone
        bool globalTrack_;

};

MuonHitCounter::MuonHitCounter(const edm::ParameterSet &iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    globalTrack_(iConfig.getParameter<bool>("useGlobalTrack"))
{

    produces<edm::ValueMap<int> >("");
    produces<edm::ValueMap<int> >("any"); 
    for(size_t j =0; j<4; ++j) {
      std::string intLabel = lexical_cast<std::string>(j+1);
      produces<edm::ValueMap<int> >("v"+intLabel);
      produces<edm::ValueMap<int> >(intLabel+"any");
      produces<edm::ValueMap<int> >("dt"+intLabel);
      produces<edm::ValueMap<int> >("dt"+intLabel+"any");
      produces<edm::ValueMap<int> >("csc"+intLabel);
      produces<edm::ValueMap<int> >("csc"+intLabel+"any");
      produces<edm::ValueMap<int> >("rpc"+intLabel);
      produces<edm::ValueMap<int> >("rpc"+intLabel+"any");
    }
}

MuonHitCounter::~MuonHitCounter() 
{
}

void
MuonHitCounter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::LogVerbatim("MuonHitCounter") <<"\n sono in MuonHitCounter !";

    edm::Handle<edm::View<reco::Muon> > src; 
    iEvent.getByLabel(src_, src);

    std::vector<int>  sumvalid(src->size(), 0);
    std::vector<int>  sumany(src->size(), 0);
    
    std::vector< std::vector<int> > valid(4,std::vector<int>(src->size(), 0));
    std::vector< std::vector<int> > any(4,std::vector<int>(src->size(), 0));
    std::vector< std::vector<int> > cscAny(4,std::vector<int>(src->size(), 0));
    std::vector< std::vector<int> > dtAny(4,std::vector<int>( src->size(), 0));
    std::vector< std::vector<int> > rpcAny(4,std::vector<int>(src->size(), 0));
    std::vector< std::vector<int> > csc(4,std::vector<int>(src->size(), 0));
    std::vector< std::vector<int> > dt(4,std::vector<int>( src->size(), 0));
    std::vector< std::vector<int> > rpc(4,std::vector<int>(src->size(), 0));
    
    for (size_t i = 0, n = src->size(); i < n; ++i) {
      const reco::Muon &mu = (*src)[i];
      reco::TrackRef track = (globalTrack_ ? mu.globalTrack() : mu.outerTrack());
      if (track.isNull()) continue;

      sumvalid[i]  = muon::muonHitCount(track,0,0,true);
      sumany[i]    = muon::muonHitCount(track,0,0,false);
      for(size_t j = 0; j < 4; ++j) {
	valid[j][i]  = muon::muonHitCount(track,0,j+1,true);    
	any[j][i]    = muon::muonHitCount(track,0,j+1,false); 
	csc[j][i]    = muon::muonHitCount(track,MuonSubdetId::CSC,j+1,true); 
	dt[j][i]     = muon::muonHitCount(track,MuonSubdetId::DT, j+1,true); 
	rpc[j][i]    = muon::muonHitCount(track,MuonSubdetId::RPC,j+1,true); 
	cscAny[j][i] = muon::muonHitCount(track,MuonSubdetId::CSC,j+1,false); 
	dtAny[j][i]  = muon::muonHitCount(track,MuonSubdetId::DT, j+1,false); 
	rpcAny[j][i] = muon::muonHitCount(track,MuonSubdetId::RPC,j+1,false); 
      }
    }
    
    writeValueMap(iEvent, src, sumvalid,   "");
    writeValueMap(iEvent, src, sumany,     "any");
    for(size_t j = 0; j<4; ++j) {
      std::string intLabel = lexical_cast<std::string>(j+1);
      writeValueMap(iEvent, src, valid[j],    "v"+intLabel);
      writeValueMap(iEvent, src, any[j],      intLabel+"any");
      writeValueMap(iEvent, src, csc[j],      "csc"+intLabel);
      writeValueMap(iEvent, src, dt[j],       "dt"+intLabel);
      writeValueMap(iEvent, src, rpc[j],      "rpc"+intLabel);
      writeValueMap(iEvent, src, cscAny[j],   "csc"+intLabel+"any");
      writeValueMap(iEvent, src, dtAny[j],    "dt"+intLabel+"any");
      writeValueMap(iEvent, src, rpcAny[j],   "rpc"+intLabel+"any");
    }
}


void
MuonHitCounter::writeValueMap(edm::Event &iEvent,
        const edm::Handle<edm::View<reco::Muon> > & handle,
        const std::vector<int> & values,
        const std::string    & label) const 
{
    using namespace edm; 
    using namespace std;

    auto_ptr<ValueMap<int> > valMap(new ValueMap<int>());
    edm::ValueMap<int>::Filler filler(*valMap);
    filler.insert(handle, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap, label);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonHitCounter);
