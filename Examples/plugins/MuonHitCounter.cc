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
// $Id: MuonHitCounter.cc,v 1.1 2010/05/13 15:31:45 gpetrucc Exp $
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
    produces<edm::ValueMap<int> >("csc"); 
    produces<edm::ValueMap<int> >("cscAny"); 
    produces<edm::ValueMap<int> >("dt"); 
    produces<edm::ValueMap<int> >("dtAny"); 
    produces<edm::ValueMap<int> >("rpc"); 
    produces<edm::ValueMap<int> >("rpcAny"); 
}

MuonHitCounter::~MuonHitCounter() 
{
}

void
MuonHitCounter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<reco::Muon> > src; 
    iEvent.getByLabel(src_, src);

 
    std::vector<int> valid(src->size(), 0);
    std::vector<int> any(src->size(), 0);
    std::vector<int> cscAny(src->size(), 0);
    std::vector<int> dtAny( src->size(), 0);
    std::vector<int> rpcAny(src->size(), 0);
    std::vector<int> csc(src->size(), 0);
    std::vector<int> dt( src->size(), 0);
    std::vector<int> rpc(src->size(), 0);
 
    for (size_t i = 0, n = src->size(); i < n; ++i) {
        const reco::Muon &mu = (*src)[i];
        reco::TrackRef track = (globalTrack_ ? mu.globalTrack() : mu.outerTrack());
        if (track.isNull()) continue;
        
        valid[i]  = muon::muonHitCount(track,0,true);    
        any[i]    = muon::muonHitCount(track,0,false); 
        csc[i]    = muon::muonHitCount(track,MuonSubdetId::CSC,true); 
        dt[i]     = muon::muonHitCount(track,MuonSubdetId::DT, true); 
        rpc[i]    = muon::muonHitCount(track,MuonSubdetId::RPC,true); 
        cscAny[i] = muon::muonHitCount(track,MuonSubdetId::CSC,false); 
        dtAny[i]  = muon::muonHitCount(track,MuonSubdetId::DT, false); 
        rpcAny[i] = muon::muonHitCount(track,MuonSubdetId::RPC,false); 
    }

    writeValueMap(iEvent, src, valid,  "");
    writeValueMap(iEvent, src, any,    "any");
    writeValueMap(iEvent, src, csc,    "csc");
    writeValueMap(iEvent, src, dt,     "dt");
    writeValueMap(iEvent, src, rpc,    "rpc");
    writeValueMap(iEvent, src, cscAny, "cscAny");
    writeValueMap(iEvent, src, dtAny,  "dtAny");
    writeValueMap(iEvent, src, rpcAny, "rpcAny");
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
