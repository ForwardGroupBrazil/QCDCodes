#ifndef MULTIJETFILTER_H
#define MULTIJETFILTER_H
#include "FWCore/Framework/interface/EDFilter.h"

template<class Jet>
class MultijetFilter : public edm::EDFilter 
{
  public:
    typedef std::vector<Jet> JetCollection; 
    explicit MultijetFilter(const edm::ParameterSet & cfg);
    virtual void beginJob();
    bool filter(edm::Event &event, const edm::EventSetup &iSetup);
    virtual void endJob();  
  private:
    edm::InputTag jets_;
    double minPt_;
    int minNjets_;
};
#include "KKousour/QCDAnalysis/plugins/MultijetFilter.icc"
#endif

