// -*- C++ -*-
//
// Package:
// Class:  
// 
/**\class

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Nov 16 16:12 (lxplus231.cern.ch)
//         Created:  Sun Nov 16 16:14:09 CET 2008
// $Id: MuonHitCounter.cc,v 1.2 2010/06/07 00:13:54 aeverett Exp $
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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"


#include "boost/lexical_cast.hpp"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using boost::lexical_cast;

//
// class decleration
class MuonCosmicMaps : public edm::EDProducer {
    public:
        explicit MuonCosmicMaps(const edm::ParameterSet&);
        ~MuonCosmicMaps();

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);

        /// Write a ValueMap<int> in the event
        template<typename T>
        void writeValueMap(edm::Event &iEvent,
                const edm::Handle<edm::View<reco::Muon> > & handle,
                const std::vector<T> & values,
                const std::string    & label) const ;

        /// The muons
        edm::InputTag src_;

        /// Use global track instead of standalone
        bool globalTrack_;
        edm::InputTag inputMuonCosmicCompatibilityValueMap_;
};

MuonCosmicMaps::MuonCosmicMaps(const edm::ParameterSet &iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    globalTrack_(iConfig.getParameter<bool>("useGlobalTrack")),
    inputMuonCosmicCompatibilityValueMap_(iConfig.getParameter<edm::InputTag>("inputMuonCosmicCompatibilityValueMap"))
{

    produces<edm::ValueMap<float> >("");
    produces<edm::ValueMap<float> >("timeCompatibility");
    produces<edm::ValueMap<float> >("backToBackCompatibility");
    produces<edm::ValueMap<float> >("overlapCompatibility");
    produces<edm::ValueMap<float> >("ipCompatibility");
    produces<edm::ValueMap<float> >("vertexCompatibility");


}

MuonCosmicMaps::~MuonCosmicMaps() 
{
}

void
MuonCosmicMaps::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::LogVerbatim("MuonCosmicMaps") <<"\n sono in MuonCosmicMaps !";

    edm::Handle<edm::View<reco::Muon> > src; 
    iEvent.getByLabel(src_, src);

    edm::Handle<edm::ValueMap<reco::MuonCosmicCompatibility> > muonCosmicCompatibilityValueMapH_;
    iEvent.getByLabel(inputMuonCosmicCompatibilityValueMap_.label(), muonCosmicCompatibilityValueMapH_);

    std::vector<float>  comp(src->size(), 0);
    std::vector<float>  time(src->size(), 0);
    std::vector<float>  btob(src->size(), 0);
    std::vector<float>  overlap(src->size(), 0);
    std::vector<float>  ip(src->size(), 0);
    std::vector<float>  vtx(src->size(), 0);

    for (size_t i = 0, n = src->size(); i < n; ++i) {
      const reco::Muon &mu = (*src)[i];
      edm::RefToBase<reco::Muon> muonRef = src->refAt(i);
      reco::TrackRef track = (globalTrack_ ? mu.globalTrack() : mu.outerTrack());
      reco::MuonCosmicCompatibility cblock = (*muonCosmicCompatibilityValueMapH_)[muonRef];

      comp[i] = cblock.cosmicCompatibility;
      time[i] = cblock.timeCompatibility;
      btob[i] = cblock.backToBackCompatibility;
      overlap[i] = cblock.overlapCompatibility;
      ip[i] = cblock.ipCompatibility;
      vtx[i] = cblock.vertexCompatibility;
    }
    
    writeValueMap(iEvent, src, comp,    "");
    writeValueMap(iEvent, src, time,    "timeCompatibility");
    writeValueMap(iEvent, src, btob,    "backToBackCompatibility");
    writeValueMap(iEvent, src, overlap, "overlapCompatibility");
    writeValueMap(iEvent, src, ip,      "ipCompatibility");
    writeValueMap(iEvent, src, vtx,     "vertexCompatibility");

}

template<typename T>
void
MuonCosmicMaps::writeValueMap(edm::Event &iEvent,
        const edm::Handle<edm::View<reco::Muon> > & handle,
        const std::vector<T> & values,
        const std::string    & label) const 
{
    using namespace edm; 
    using namespace std;

    auto_ptr<ValueMap<T> > valMap(new ValueMap<T>());
    typename edm::ValueMap<T>::Filler filler(*valMap);
    filler.insert(handle, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap, label);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonCosmicMaps);
